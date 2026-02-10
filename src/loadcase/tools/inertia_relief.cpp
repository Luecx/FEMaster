//
// Created by f_eggers on 10.02.2026.
//

#include "inertia_relief.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"
#include "../../model/element/element_structural.h"
#include "../../bc/load.h"

#include <Eigen/Cholesky>
#include <iomanip>
#include <utility>
#include <algorithm>

namespace fem {

void apply_inertia_relief(model::ModelData& model_data, model::Field& global_load_mat) {
    constexpr Precision kIreg   = static_cast<Precision>(1e-12);
    constexpr Precision kRelTol = static_cast<Precision>(1e-8);
    constexpr Precision kAbsTol = static_cast<Precision>(1e-10);

    using std::setw;
    using std::setprecision;

    constexpr int    w = 12;
    constexpr int    p = 3;

    auto& data = model_data;

    logging::error(data.positions != nullptr, "InertiaRelief: positions not initialized");
    auto& pos = *data.positions;

    logging::error(global_load_mat.rows == pos.rows,
                   "InertiaRelief: global_load_mat.rows != positions.rows (",
                   global_load_mat.rows, " vs ", pos.rows, ")");

    // ---------------------------------------------------------------------
    // Helper: compute resultant force & moment about center c
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
    // Mass and center of gravity
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

        m     += se->integrate_scalar_field(true, [](const Vec3&) { return Precision(1); });
        c_num += se->integrate_vector_field(true, [](const Vec3& x) { return x; });
    }

    logging::error(m > Precision(0), "InertiaRelief: total mass is zero -> cannot apply inertia relief");
    const Vec3 c = c_num / m;

    // ---------------------------------------------------------------------
    // Initial resultants
    // ---------------------------------------------------------------------
    const auto [Fext, Mext] = compute_resultants(global_load_mat, c);

    // ---------------------------------------------------------------------
    // Inertia tensor about CoG
    // ---------------------------------------------------------------------
    Mat3 I = Mat3::Zero();

    for (auto& eptr : data.elements) {
        if (!eptr) {
            continue;
        }
        auto* se = eptr->as<model::StructuralElement>();
        if (!se) {
            continue;
        }

        I += se->integrate_tensor_field(true, [c](const Vec3& x) {
            const Vec3       r = x - c;
            const Precision r2 = r.squaredNorm();
            const Mat3 integrand =
                r2 * Mat3::Identity() - (r * r.transpose());
            return integrand;
        });
    }

    // ---------------------------------------------------------------------
    // Solve rigid-body accelerations
    // (sign convention consistent with existing InertialLoad)
    // ---------------------------------------------------------------------
    const Vec3 a0 = Fext / m;

    Vec3 alpha = Vec3::Zero();
    {
        Eigen::SelfAdjointEigenSolver<Mat3> es(I);
        logging::error(es.info() == Eigen::Success,
                       "InertiaRelief: eigen decomposition of inertia tensor failed");

        const Vec3 evals = es.eigenvalues();
        const Mat3 evecs = es.eigenvectors();

        const Precision eval_max = evals.cwiseAbs().maxCoeff();
        if (!(std::isfinite(eval_max) && eval_max > Precision(0))) {
            logging::info(true, I);
            logging::error(std::isfinite(eval_max) && eval_max > Precision(0),
                           "InertiaRelief: inertia tensor is zero/invalid");
        }

        // Cutoff: treat very small principal inertias as unconstrained rotations
        // (i.e. alpha component in that direction is set to 0).
        const Precision eps = std::max(kIreg, Precision(1e-12) * eval_max);

        Vec3 inv = Vec3::Zero();
        for (int i = 0; i < 3; ++i) {
            if (std::abs(evals(i)) > eps) {
                inv(i) = Precision(1) / evals(i);
            } else {
                inv(i) = Precision(0);
            }
        }

        alpha = evecs * inv.asDiagonal() * evecs.transpose() * Mext;
    }

    logging::error(alpha.allFinite() && a0.allFinite(),
                   "InertiaRelief: computed accelerations contain NaN/Inf");

    // ---------------------------------------------------------------------
    // Logging: 3 rows, 6 columns
    // ---------------------------------------------------------------------
    logging::info(true, "InertiaRelief:");
    logging::info(true, "              ",
                  std::setw(w), "x" , std::setw(w), "y" , std::setw(w), "z" ,
                  std::setw(w), "rx", std::setw(w), "ry", std::setw(w), "rz");

    logging::info(true, "initial F/M  :", setprecision(p),
                  std::setw(w), Fext(0), std::setw(w), Fext(1), std::setw(w), Fext(2),
                  std::setw(w), Mext(0), std::setw(w), Mext(1), std::setw(w), Mext(2));

    logging::info(true, "using a/α    :", setprecision(p),
                  std::setw(w), a0(0), std::setw(w), a0(1), std::setw(w), a0(2),
                  std::setw(w), alpha(0), std::setw(w), alpha(1), std::setw(w), alpha(2));

    // ---------------------------------------------------------------------
    // Apply inertial load directly to a temporary field and subtract it
    // ---------------------------------------------------------------------
    model::Field inertial_mat = data.create_field_("_TEMP", model::FieldDomain::NODE, 6);
    inertial_mat.set_zero();

    fem::bc::InertialLoad ir;
    ir.center     = c;
    ir.center_acc = a0;
    ir.alpha      = alpha;
    ir.omega      = Vec3::Zero();
    ir.region     = model_data.elem_sets.all();

    ir.apply(data, inertial_mat, Precision(0));

    global_load_mat += inertial_mat;

    // ---------------------------------------------------------------------
    // Resulting residuals and sanity check
    // ---------------------------------------------------------------------
    const auto [Fsum, Msum] = compute_resultants(global_load_mat, c);

    logging::info(true, setprecision(p),
              "resulting F/M:",
              setw(w), Fsum(0), setw(w), Fsum(1), setw(w), Fsum(2),
              setw(w), Msum(0), setw(w), Msum(1), setw(w), Msum(2));

    const Precision tolF = std::max(kAbsTol, kRelTol * (Precision(1) + Fext.cwiseAbs().maxCoeff()));
    const Precision tolM = std::max(kAbsTol, kRelTol * (Precision(1) + Mext.cwiseAbs().maxCoeff()));

    logging::error(Fsum.template lpNorm<Eigen::Infinity>() <= tolF &&
                   Msum.template lpNorm<Eigen::Infinity>() <= tolM,
                   "InertiaRelief: residual balance too large (|F|=", Fsum.norm(), ", |M|=", Msum.norm(), ")");
}

} // namespace fem