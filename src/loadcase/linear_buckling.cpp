// loadcases/linear_buckling.cpp
#include "linear_buckling.h"

#include <algorithm>  // std::sort
#include <iomanip>
#include <limits>

#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/extract_scaled_row_sum.h"
#include "../solve/eigen.h"  // new API: eigs(...) -> std::vector<EigenValueVectorPair>
#include "../reader/write_mtx.h"

namespace fem { namespace loadcase {

struct BucklingMode {
    Precision     lambda;            // buckling factor
    DynamicVector mode_shape;        // reduced (includes Lagrange rows at tail)
    DynamicMatrix mode_shape_mat;    // expanded to node x dof layout

    explicit BucklingMode(Precision lam, DynamicVector v)
        : lambda(lam), mode_shape(std::move(v)) {}

    void expand_mode_shape(int lagrange_dofs,
                           const DynamicVector& active_lhs_vec,
                           const IndexMatrix&   active_dof_idx_mat)
    {
        DynamicVector mode_shape_active = mode_shape.head(mode_shape.size() - lagrange_dofs);
        auto expanded_mode_shape        = mattools::expand_vec_to_vec(mode_shape_active, active_lhs_vec);
        mode_shape_mat                  = mattools::expand_vec_to_mat(active_dof_idx_mat, expanded_mode_shape);
    }

    bool operator<(const BucklingMode& o) const { return lambda < o.lambda; }
};

// Convert solver eigenpairs to sorted buckling modes
static std::vector<BucklingMode>
    make_and_sort_modes(const std::vector<solver::EigenValueVectorPair>& pairs)
{
    std::vector<BucklingMode> out;
    out.reserve(pairs.size());
    for (const auto& p : pairs)
        out.emplace_back(p.value, p.vector);
    std::sort(out.begin(), out.end());
    return out;
}

/// Check symmetry and compute eigenvalues of a (small/medium) sparse matrix.
/// Only use this as a debug check (will dense-convert internally).
static void check_matrix_spectrum(const Eigen::SparseMatrix<double>& A, int num_lowest = 10)
{
    Eigen::MatrixXd Ad = Eigen::MatrixXd(A);

    // Symmetry check
    Eigen::MatrixXd diff = Ad - Ad.transpose();
    double sym_norm = diff.norm();
    double a_norm   = Ad.norm();
    double rel_sym  = (a_norm > 0) ? sym_norm / a_norm : sym_norm;

    std::cout << "\n[DEBUG] Matrix spectrum check\n";
    std::cout << "  Size           : " << Ad.rows() << " x " << Ad.cols() << "\n";
    std::cout << "  Symmetry error : " << sym_norm
              << " (relative " << rel_sym << ")\n";

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ad);
    if (es.info() != Eigen::Success) {
        std::cout << "  [ERROR] Eigen decomposition failed.\n";
        return;
    }

    Eigen::VectorXd eigs = es.eigenvalues();
    std::cout << "  Lowest " << num_lowest << " eigenvalues:\n";
    int n = std::min<int>(num_lowest, eigs.size());
    for (int i = 0; i < n; ++i)
        std::cout << "    λ[" << i << "] = " << eigs[i] << "\n";
    std::cout.flush();
}

// --- Simple, compact on-screen summary like eigenfrequency (Idx, λ only) ---
static void display_buckling_results(const std::vector<BucklingMode>& modes,
                                     double sigma_used,
                                     int    ncv_used,
                                     int    k_requested)
{
    logging::info(true, "");
    logging::info(true, "Buckling summary");
    logging::up();
    logging::info(true, "requested modes k : ", k_requested);
    logging::info(true, "ncv               : ", ncv_used);
    logging::info(true, "shift σ           : ", std::scientific, std::setprecision(3), sigma_used);
    logging::down();

    logging::info(true, std::setw(6), "Idx", std::setw(24), "Buckling factor λ");
    for (size_t i = 0; i < modes.size(); ++i) {
        logging::info(true,
                      std::setw(6), i + 1,
                      std::setw(24), std::fixed, std::setprecision(9), Precision(modes[i].lambda));
    }
}

LinearBuckling::LinearBuckling(ID id, reader::Writer* writer, model::Model* model, int numEigenvalues)
    : LoadCase(id, writer, model), num_eigenvalues(numEigenvalues) {}

void LinearBuckling::run() {
    logging::info(true, "\n\n================================================================================");
    logging::info(true, "LINEAR BUCKLING");
    logging::info(true, "================================================================================\n");

    // 0) Sections/materials
    m_model->assign_sections();

    // 1) DOF indexing + supports + loads
    auto active_dof_idx_mat = Timer::measure(
        [&]{ return m_model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix" );

    auto [global_supp_mat, global_supp_eqs] = Timer::measure(
        [&]{ return m_model->build_support_matrix(supps); },
        "building global support matrix" );

    auto global_load_mat = Timer::measure(
        [&]{ return m_model->build_load_matrix(loads); },
        "building global load matrix" );

    // 2) Assemble linear K and constraints (like LinearStatic)
    auto K_act = Timer::measure(
        [&]{ return m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing active stiffness matrix" );

    Precision kchar = K_act.diagonal().mean();

    auto C_act = Timer::measure(
        [&]{ return m_model->build_constraint_matrix(active_dof_idx_mat, global_supp_eqs, kchar); },
        "constructing active Lagrangian matrix" );

    const int m = K_act.rows();     // active DOFs
    const int n = C_act.rows();     // Lagrange rows

    // 3) Assemble augmented K̂ = [K  Cᵀ; C  0]
    auto K_hat = Timer::measure(
        [&]() {
            SparseMatrix A(m + n, m + n);
            TripletList T; T.reserve(K_act.nonZeros() + 2*C_act.nonZeros() + n);

            // K
            for (int k = 0; k < K_act.outerSize(); ++k)
                for (SparseMatrix::InnerIterator it(K_act, k); it; ++it)
                    T.emplace_back(it.row(), it.col(), it.value());

            // C blocks
            for (int k = 0; k < C_act.outerSize(); ++k)
                for (SparseMatrix::InnerIterator it(C_act, k); it; ++it) {
                    T.emplace_back(m + it.row(), it.col(), it.value());
                    T.emplace_back(it.col(), m + it.row(), it.value());
                }

            // small regularization in bottom-right (Lagrange block)
            for (int i = 0; i < n; ++i)
                T.emplace_back(m + i, m + i, -kchar / 1e6);

            A.setFromTriplets(T.begin(), T.end());
            return A;
        },
        "assembling full lhs matrix K_hat" );

    // 4) Reduced RHS/LHS for the static pre-stress solve (no output)
    auto active_rhs_vec = Timer::measure(
        [&]{ return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load vector (RHS)" );

    auto active_lhs_vec = Timer::measure(
        [&]{ return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_supp_mat); },
        "reducing support vector (LHS)" );

    DynamicVector full_rhs_vec(m + n);
    DynamicVector full_lhs_vec(m + n);
    full_rhs_vec << active_rhs_vec, DynamicVector::Zero(n);
    full_lhs_vec << active_lhs_vec, DynamicVector::Constant(n, std::numeric_limits<Precision>::quiet_NaN());

    // Implicit RHS (like LinearStatic)
    auto implicit_rhs = Timer::measure(
        [&]{ return mattools::extract_scaled_row_sum(K_hat, full_lhs_vec); },
        "computing implicit load vector" );
    full_rhs_vec += implicit_rhs;

    // Reduce K̂ and RHS
    auto K_hat_red = Timer::measure(
        [&]{ return mattools::reduce_mat_to_mat(K_hat, full_lhs_vec); },
        "reducing K_hat to solver-ready form" );

    auto rhs_red = Timer::measure(
        [&]{ return mattools::reduce_vec_to_vec(full_rhs_vec, full_lhs_vec); },
        "reducing load vector to solver-ready form" );

    K_hat_red.makeCompressed();

    // 5) Solve static displacements (only to form IP stresses)
    auto sol = Timer::measure(
        [&]{ return solve(device, method, K_hat_red, rhs_red); },
        "solving linear system for pre-stress" );

    // Expand back to full DOF vector [u; λ]
    DynamicVector u_active = sol.head(sol.rows() - n);
    DynamicVector u_full   = mattools::expand_vec_to_vec(u_active, active_lhs_vec);
    DynamicMatrix U_mat    = mattools::expand_vec_to_mat(active_dof_idx_mat, u_full);

    // 6) Compute IP stresses for geometric stiffness
    IPData ip_stress;
    {
        NodeData U_row = U_mat;
        IPData ip_strain_unused;
        std::tie(ip_stress, ip_strain_unused) =
            Timer::measure([&]{ return m_model->compute_ip_stress_strain(U_row); },
                           "computing IP stress/strain for Kg");
    }

    // 7) Assemble geometric stiffness Kg (active) and augment to K̂g like K̂
    auto Kg_act = Timer::measure(
        [&]{ return m_model->build_geom_stiffness_matrix(active_dof_idx_mat, ip_stress); },
        "assembling active geometric stiffness matrix" );

    auto Kg_hat = Timer::measure(
        [&]() {
            SparseMatrix G(m + n, m + n);
            TripletList T; T.reserve(Kg_act.nonZeros());
            // put Kg in top-left (Λ blocks remain zero)
            for (int k = 0; k < Kg_act.outerSize(); ++k)
                for (SparseMatrix::InnerIterator it(Kg_act, k); it; ++it)
                    T.emplace_back(it.row(), it.col(), it.value());
            G.setFromTriplets(T.begin(), T.end());
            return G;
        },
        "assembling Kghat" );

    // Reduce Kg_hat consistent with constraints
    auto Kg_hat_red = Timer::measure(
        [&]{ return mattools::reduce_mat_to_mat(Kg_hat, full_lhs_vec); },
        "reducing Kghat to solver-ready form" );

    K_hat_red.makeCompressed();
    Kg_hat_red.makeCompressed();

    // Overview
    logging::info(true, "");
    logging::info(true, "Overview");
    logging::up();
    logging::info(true, "system total DOFs : ", active_dof_idx_mat.maxCoeff() + 1);
    logging::info(true, "lagrange DOFs     : ", n);
    logging::info(true, "final DOFs        : ", rhs_red.rows());
    logging::info(true, "active system DOFs: ", rhs_red.rows() - n);
    logging::down();

    // 8) Buckling EVP
    // Classical formulation: K φ + λ Kg φ = 0  ⇔  K φ = (-λ) Kg φ.
    // If Kg is assembled in the classical sign convention, use B = -Kg_hat to obtain positive buckling factors.
    SparseMatrix B = -Kg_hat_red;

    const Index active_system_dofs = rhs_red.rows() - n;
    const Index neigs = std::min<Index>(num_eigenvalues, std::max<Index>(1, active_system_dofs));

    // --- Solver options (choose σ as you like; 1e-9 recommended by you) ---
    solver::EigenOpts opts;
    opts.mode = solver::EigenMode::Buckling;
    opts.sigma = 1e-3; // <<< your tiny shift; change here if desired
    opts.sort = solver::EigenOpts::Sort::SmallestAlge; // smallest λ first
    const int ncv = std::min<int>(static_cast<int>(K_hat_red.rows()),
                                  std::max<int>(static_cast<int>(3 * neigs + 20), static_cast<int>(neigs + 2)));
    opts.ncv = ncv;

    // Optional: debug spectra of K and B
    // check_matrix_spectrum(K_hat_red);
    // check_matrix_spectrum(B);

    // Solve with Spectra Buckling mode via new API (back-transforms to original λ)
    auto eig_pairs = Timer::measure(
        [&]{
            return solver::eigs(solver::CPU, K_hat_red, B,
                                static_cast<int>(neigs), opts);
        },
        "solving generalized EVP for buckling"
    );

    // 9) Pair/sort and expand
    auto modes = make_and_sort_modes(eig_pairs);
    for (auto& m : modes)
        m.expand_mode_shape(n, active_lhs_vec, active_dof_idx_mat);

    // --- NEW: compact on-screen summary ---
    display_buckling_results(modes, opts.sigma, ncv, static_cast<int>(neigs));

    // 10) Write results: buckling factors + mode shapes
    m_writer->add_loadcase(m_id);
    {
        DynamicVector lambdas(modes.size());
        for (size_t i = 0; i < modes.size(); ++i) {
            lambdas(i) = modes[i].lambda;
            m_writer->write_eigen_matrix(modes[i].mode_shape_mat, "BUCKLING_MODE_" + std::to_string(i+1));
        }
        m_writer->write_eigen_matrix(DynamicMatrix(lambdas), "BUCKLING_FACTORS");
    }

    logging::info(true, "Buckling analysis completed.");
}

}} // namespace fem::loadcase
