// loadcases/linear_buckling.cpp
#include "linear_buckling.h"

#include <algorithm>  // std::sort
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/extract_scaled_row_sum.h"
#include "../solve/eigen.h"

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

static std::vector<BucklingMode>
pair_and_sort(const DynamicVector& lambdas, const DynamicMatrix& vecs)
{
    std::vector<BucklingMode> out;
    out.reserve(lambdas.size());
    for (int i = 0; i < lambdas.size(); ++i)
        out.emplace_back(lambdas[i], vecs.col(i));
    std::sort(out.begin(), out.end());
    return out;
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

    // 6) Compute IP stresses for geometric stiffness (no need to pass an enumeration)
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

            // put Kg in top-left
            for (int k = 0; k < Kg_act.outerSize(); ++k)
                for (SparseMatrix::InnerIterator it(Kg_act, k); it; ++it)
                    T.emplace_back(it.row(), it.col(), it.value());

            // Lagrange blocks remain zero
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
    logging::info(true, "\nOverview");
    logging::up();
    logging::info(true, "system total DOFs : ", active_dof_idx_mat.maxCoeff() + 1);
    logging::info(true, "lagrange DOFs     : ", n);
    logging::info(true, "final DOFs        : ", rhs_red.rows());
    logging::down();

    // 8) Buckling EVP:  K_hat * φ = λ * (-Kg_hat) * φ  => use B = (-Kg_hat)
    const Index active_system_dofs = rhs_red.rows() - n;
    const Index neigs = std::min<Index>(num_eigenvalues, std::max<Index>(1, active_system_dofs));

    auto sym = [](SparseMatrix& A){
        A.makeCompressed();
    };
    sym(K_hat_red);
    sym(Kg_hat_red);

    SparseMatrix B = Kg_hat_red;

    auto eig = Timer::measure(
        [&]{ return solver::compute_eigenvalues(solver::CPU, K_hat_red, B, neigs, true); },
        "solving generalized EVP for buckling" );

    // 9) Pair/sort and expand
    auto modes = pair_and_sort(eig.first, eig.second);
    for (auto& m : modes)
        m.expand_mode_shape(n, active_lhs_vec, active_dof_idx_mat);

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
