//
// Created by f_eggers on 10.02.2026.
//

#include "rebalance_loads.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"
#include "../../model/element/element_structural.h"

#include <Eigen/Dense>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cmath>

namespace fem {

void rebalance_loads(model::ModelData& model_data, model::Field& global_load_mat) {
    constexpr Precision kRelTol = static_cast<Precision>(1e-10);
    constexpr Precision kAbsTol = static_cast<Precision>(1e-12);
    constexpr Precision kDlsReg = static_cast<Precision>(1e-14);

    using std::setw;
    using std::setprecision;

    constexpr int w = 12;
    constexpr int p = 3;

    auto& data = model_data;

    logging::error(data.positions != nullptr, "RigidBodyRebalancing: positions not initialized");
    auto& pos = *data.positions;

    logging::error(global_load_mat.rows == pos.rows,
                   "RigidBodyRebalancing: global_load_mat.rows != positions.rows (",
                   global_load_mat.rows, " vs ", pos.rows, ")");

    logging::error(global_load_mat.components >= 3,
                   "RigidBodyRebalancing: global_load_mat.components must be >= 3");

    // ---------------------------------------------------------------------
    // Helper: compute resultant force & moment about reference point c
    // Includes nodal moments if present (components >= 6).
    // ---------------------------------------------------------------------
    const auto compute_resultants = [&](const model::Field& L, const Vec3& c) -> std::pair<Vec3, Vec3> {
        Vec3 F = Vec3::Zero();
        Vec3 M = Vec3::Zero();

        const bool has_moments = (L.components >= 6);
        for (Index i = 0; i < L.rows; ++i) {
            const Vec3 fi{L(i, 0), L(i, 1), L(i, 2)};
            Vec3 mi = Vec3::Zero();
            if (has_moments) {
                mi = Vec3{L(i, 3), L(i, 4), L(i, 5)};
            }

            const Vec3 xi = pos.row_vec3(i);
            F += fi;
            M += mi + (xi - c).cross(fi);
        }
        return {F, M};
    };

    // ---------------------------------------------------------------------
    // Mass and center of gravity (re-using your structural element integration)
    // ---------------------------------------------------------------------
    Precision m = Precision(0);
    Vec3 c_num  = Vec3::Zero();

    for (auto& eptr : data.elements) {
        if (!eptr) {
            continue;
        }
        auto* se = eptr->as<model::StructuralElement>();
        if (!se) {
            continue;
        }

        // NOTE: your integrate_scalar_field(true, ...) seems to integrate rho*dV
        // or dV depending on implementation; inertia_relief uses the same.
        m     += se->integrate_scalar_field(true, [](const Vec3&) { return Precision(1); });
        c_num += se->integrate_vector_field(true, [](const Vec3& x) { return x; });
    }

    logging::error(m > Precision(0),
                   "RigidBodyRebalancing: total mass is zero -> cannot determine CoG");
    const Vec3 c = c_num / m;

    // ---------------------------------------------------------------------
    // Compute initial residual resultants
    // ---------------------------------------------------------------------
    const auto [F0, M0] = compute_resultants(global_load_mat, c);

    const Precision tolF = std::max(kAbsTol, kRelTol * (Precision(1) + F0.cwiseAbs().maxCoeff()));
    const Precision tolM = std::max(kAbsTol, kRelTol * (Precision(1) + M0.cwiseAbs().maxCoeff()));

    // ---------------------------------------------------------------------
    // Logging header (same shape as InertiaRelief)
    // ---------------------------------------------------------------------
    logging::info(true, "RigidBodyRebalancing:");
    logging::info(true, "              ",
                  std::setw(w), "x" , std::setw(w), "y" , std::setw(w), "z" ,
                  std::setw(w), "rx", std::setw(w), "ry", std::setw(w), "rz");

    logging::info(true, "initial F/M  :", setprecision(p),
                  std::setw(w), F0(0), std::setw(w), F0(1), std::setw(w), F0(2),
                  std::setw(w), M0(0), std::setw(w), M0(1), std::setw(w), M0(2));

    // Early out if already balanced
    if (F0.template lpNorm<Eigen::Infinity>() <= tolF &&
        M0.template lpNorm<Eigen::Infinity>() <= tolM) {
        logging::info(true, "RigidBodyRebalancing: already balanced within tolerance.");
        return;
    }

    // ---------------------------------------------------------------------
    // Build a 6D rigid-load basis at nodal level (forces only):
    //   T_x, T_y, T_z: uniform force directions
    //   R_x, R_y, R_z: force field that produces a pure moment about c
    // ---------------------------------------------------------------------
    const Vec3 ex{Precision(1), Precision(0), Precision(0)};
    const Vec3 ey{Precision(0), Precision(1), Precision(0)};
    const Vec3 ez{Precision(0), Precision(0), Precision(1)};

    // TODO: consider nodal lumped masses for wgt(i) to get "physically nicer" corrections.
    Eigen::VectorXd wgt(pos.rows);
    wgt.setOnes();

    // Helper: compute resultants of a basis load j (j=0..5) without allocating a Field.
    const auto basis_resultants = [&](int j) -> std::pair<Vec3, Vec3> {
        Vec3 F = Vec3::Zero();
        Vec3 M = Vec3::Zero();

        const Vec3 axis = (j == 0) ? ex : (j == 1) ? ey : (j == 2) ? ez
                         : (j == 3) ? ex : (j == 4) ? ey : ez;

        for (Index i = 0; i < pos.rows; ++i) {
            const Vec3 xi = pos.row_vec3(i);
            const Vec3 r  = xi - c;

            Vec3 fi = Vec3::Zero();
            if (j < 3) {
                // Translation basis: fi = w_i * e_k
                fi = wgt(i) * axis;
            } else {
                // Rotation basis: fi = w_i * (axis x r)
                fi = wgt(i) * axis.cross(r);
            }

            F += fi;
            M += r.cross(fi);
        }

        return {F, M};
    };

    // Assemble A (6x6) mapping basis coefficients -> (F, M)
    Eigen::Matrix<Precision, 6, 6> A;
    A.setZero();
    for (int j = 0; j < 6; ++j) {
        const auto [Fj, Mj] = basis_resultants(j);
        A.template block<3,1>(0, j) = Fj;
        A.template block<3,1>(3, j) = Mj;
    }

    // RHS: cancel initial residuals
    Eigen::Matrix<Precision, 6, 1> g;
    g.template segment<3>(0) = F0;
    g.template segment<3>(3) = M0;

    // Solve A c = -g (robust: try direct, fallback to damped LS)
    Eigen::Matrix<Precision, 6, 1> coeff = Eigen::Matrix<Precision, 6, 1>::Zero();

    bool solved = false;
    {
        Eigen::FullPivLU<Eigen::Matrix<Precision, 6, 6>> lu(A);
        if (lu.isInvertible()) {
            coeff  = lu.solve(-g);
            solved = coeff.allFinite();
        }
    }

    if (!solved) {
        // Damped least squares: (A^T A + λ I) c = -A^T g
        const Eigen::Matrix<Precision, 6, 6> AtA = A.transpose() * A;

        Precision diag_max = Precision(0);
        for (int i = 0; i < 6; ++i) {
            diag_max = std::max(diag_max, std::abs(AtA(i, i)));
        }
        const Precision lambda = std::max(kDlsReg, Precision(1e-12) * (Precision(1) + diag_max));

        Eigen::Matrix<Precision, 6, 6> H = AtA + lambda * Eigen::Matrix<Precision, 6, 6>::Identity();
        Eigen::Matrix<Precision, 6, 1> rhs = -A.transpose() * g;

        Eigen::LDLT<Eigen::Matrix<Precision, 6, 6>> ldlt(H);
        logging::error(ldlt.info() == Eigen::Success,
                       "RigidBodyRebalancing: DLS LDLT failed");
        coeff = ldlt.solve(rhs);

        logging::error(coeff.allFinite(),
                       "RigidBodyRebalancing: DLS coefficients contain NaN/Inf");
    }

    // ---------------------------------------------------------------------
    // Compute (ΔF, ΔM) induced by the chosen coefficients (purely for logging)
    // ---------------------------------------------------------------------
    Eigen::Matrix<Precision, 6, 1> d = A * coeff; // should be approx -g
    const Vec3 dF = d.template segment<3>(0);
    const Vec3 dM = d.template segment<3>(3);

    logging::info(true, "coeff (Tx,Ty,Tz,Rx,Ry,Rz):", setprecision(p),
                  setw(w), coeff(0), setw(w), coeff(1), setw(w), coeff(2),
                  setw(w), coeff(3), setw(w), coeff(4), setw(w), coeff(5));

    logging::info(true, "applied ΔF/M :", setprecision(p),
                  std::setw(w), dF(0), std::setw(w), dF(1), std::setw(w), dF(2),
                  std::setw(w), dM(0), std::setw(w), dM(1), std::setw(w), dM(2));

    // ---------------------------------------------------------------------
    // Apply correction to global_load_mat (forces only)
    // ---------------------------------------------------------------------
    for (Index i = 0; i < pos.rows; ++i) {
        const Vec3 xi = pos.row_vec3(i);
        const Vec3 r  = xi - c;

        Vec3 df = Vec3::Zero();

        // Translation components
        df += wgt(i) * coeff(0) * ex;
        df += wgt(i) * coeff(1) * ey;
        df += wgt(i) * coeff(2) * ez;

        // Rotation components
        df += wgt(i) * coeff(3) * ex.cross(r);
        df += wgt(i) * coeff(4) * ey.cross(r);
        df += wgt(i) * coeff(5) * ez.cross(r);

        global_load_mat(i, 0) += df(0);
        global_load_mat(i, 1) += df(1);
        global_load_mat(i, 2) += df(2);
    }

    // ---------------------------------------------------------------------
    // Post-check
    // ---------------------------------------------------------------------
    const auto [F1, M1] = compute_resultants(global_load_mat, c);

    logging::info(true, setprecision(p),
                  "resulting F/M:", setw(w), F1(0), setw(w), F1(1), setw(w), F1(2),
                  setw(w), M1(0), setw(w), M1(1), setw(w), M1(2));

    const Precision tolF_post = std::max(kAbsTol, kRelTol * (Precision(1) + F0.cwiseAbs().maxCoeff()));
    const Precision tolM_post = std::max(kAbsTol, kRelTol * (Precision(1) + M0.cwiseAbs().maxCoeff()));

    logging::error(F1.template lpNorm<Eigen::Infinity>() <= tolF_post &&
                   M1.template lpNorm<Eigen::Infinity>() <= tolM_post,
                   "RigidBodyRebalancing: residual balance too large (|F|=", F1.norm(), ", |M|=", M1.norm(), ")");
}

} // namespace fem