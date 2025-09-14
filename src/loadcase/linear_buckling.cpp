/******************************************************************************
* @file linear_buckling.cpp
* @brief Implementation of the LinearBuckling load case for performing linear
* buckling analysis.
*
* @details
* This file implements the linear buckling analysis based on the finite element
* method. The procedure involves:
* 1. Building the stiffness matrix with constraints.
* 2. Solving a static preload step to obtain stresses.
* 3. Assembling the geometric stiffness matrix from stresses.
* 4. Solving the generalized eigenvalue problem:
*
*      K φ = λ (−Kg) φ
*
*    where K is the structural stiffness matrix and Kg is the geometric stiffness
*    matrix. The resulting eigenvalues λ are the buckling factors, and the eigenvectors
*    φ are the corresponding buckling mode shapes.
*
* A coarse eigenvalue estimate is computed cheaply (via a Rayleigh quotient using
* the prebuckling displacement) to set the shift σ for the eigenvalue solver.
*
* @date Created on 12.09.2025
******************************************************************************/

#include "linear_buckling.h"

#include <algorithm>  // std::sort
#include <iomanip>
#include <limits>

#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/extract_scaled_row_sum.h"
#include "../solve/eigen.h"   // eigs(...) -> std::vector<EigenValueVectorPair>
#include "../reader/write_mtx.h"

namespace fem { namespace loadcase {

//==============================================================================
// BucklingMode: holds a single eigenpair (buckling factor + mode shape)
//==============================================================================
struct BucklingMode {
   Precision     lambda;            ///< Buckling factor
   DynamicVector mode_shape;        ///< Reduced mode shape (includes Lagrange DOFs at tail)
   DynamicMatrix mode_shape_mat;    ///< Expanded to node × DOF layout

   explicit BucklingMode(Precision lam, DynamicVector v)
       : lambda(lam), mode_shape(std::move(v)) {}

   /// Expand reduced mode shape to full node × DOF layout
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

//==============================================================================
// Helper: convert raw solver eigenpairs into sorted BucklingModes
//==============================================================================
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

//==============================================================================
// Helper: print compact buckling summary
//==============================================================================
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
                     std::setw(24), std::fixed, std::setprecision(9),
                     Precision(modes[i].lambda));
   }
}

//==============================================================================
// Helper: estimate the first eigenvalue with a cheap Rayleigh quotient
//==============================================================================
static double estimate_first_eigenvalue(const SparseMatrix& K,
                                       const SparseMatrix& Kg,
                                       const DynamicVector& trial_u)
{
   double lambda_est = 0.0;

   DynamicVector Ku = K * trial_u;
   DynamicVector Bu = (-Kg) * trial_u;   // use −Kg to ensure positive denominator

   double num = trial_u.dot(Ku);
   double den = trial_u.dot(Bu);

   if (den > 1e-14) {
       lambda_est = num / den;
   } else {
       // fallback: diagonal ratio
       double min_ratio = std::numeric_limits<double>::infinity();
       for (int k = 0; k < K.rows(); ++k) {
           double bii = -Kg.coeff(k, k);
           if (bii > 1e-14)
               min_ratio = std::min(min_ratio, K.coeff(k, k) / bii);
       }
       lambda_est = (min_ratio < std::numeric_limits<double>::infinity())
                        ? min_ratio : 100.0; // fallback constant
   }
   return lambda_est;
}

//==============================================================================
// LinearBuckling class
//==============================================================================

LinearBuckling::LinearBuckling(ID id, reader::Writer* writer, model::Model* model, int numEigenvalues)
   : LoadCase(id, writer, model), num_eigenvalues(numEigenvalues) {}

/**
* @brief Executes the linear buckling analysis.
*
* Steps:
* 1. Assign sections and materials.
* 2. Build DOF indexing, supports, and loads.
* 3. Assemble the structural stiffness matrix and constraint matrix.
* 4. Assemble augmented system with Lagrange multipliers.
* 5. Solve static preload (linear system) to obtain displacements.
* 6. Compute stresses at integration points.
* 7. Assemble geometric stiffness matrix.
* 8. Estimate smallest eigenvalue to choose shift σ.
* 9. Solve the generalized eigenvalue problem for buckling.
* 10. Expand mode shapes and write results.
*/
void LinearBuckling::run() {
   // =========================================================================
   // Begin logging
   // =========================================================================
   logging::info(true, "");
   logging::info(true, "================================================================================");
   logging::info(true, "LINEAR BUCKLING ANALYSIS");
   logging::info(true, "================================================================================");
   logging::info(true, "");

   // -------------------------------------------------------------------------
   // 0) Sections and materials
   // -------------------------------------------------------------------------
   m_model->assign_sections();

   // -------------------------------------------------------------------------
   // 1) DOF indexing + supports + loads
   // -------------------------------------------------------------------------
   auto active_dof_idx_mat = Timer::measure(
       [&]{ return m_model->build_unconstrained_index_matrix(); },
       "generating active_dof_idx_mat index matrix" );

   auto [global_supp_mat, global_supp_eqs] = Timer::measure(
       [&]{ return m_model->build_support_matrix(supps); },
       "building global support matrix" );

   auto global_load_mat = Timer::measure(
       [&]{ return m_model->build_load_matrix(loads); },
       "building global load matrix" );

   // -------------------------------------------------------------------------
   // 2) Assemble stiffness and constraints
   // -------------------------------------------------------------------------
   auto K_act = Timer::measure(
       [&]{ return m_model->build_stiffness_matrix(active_dof_idx_mat); },
       "constructing active stiffness matrix" );

   Precision kchar = K_act.diagonal().mean();

   auto C_act = Timer::measure(
       [&]{ return m_model->build_constraint_matrix(active_dof_idx_mat, global_supp_eqs, kchar); },
       "constructing active Lagrangian matrix" );

   const int m = K_act.rows();     // active DOFs
   const int n = C_act.rows();     // Lagrange rows

   // -------------------------------------------------------------------------
   // 3) Assemble augmented K̂ = [K  Cᵀ; C  0]
   // -------------------------------------------------------------------------
   auto K_hat = Timer::measure(
       [&]() {
           SparseMatrix A(m + n, m + n);
           TripletList T; T.reserve(K_act.nonZeros() + 2*C_act.nonZeros() + n);

           // K block
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

   // -------------------------------------------------------------------------
   // 4) RHS/LHS vectors for static solve
   // -------------------------------------------------------------------------
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

   // Implicit RHS (constraint enforcement)
   auto implicit_rhs = Timer::measure(
       [&]{ return mattools::extract_scaled_row_sum(K_hat, full_lhs_vec); },
       "computing implicit load vector" );
   full_rhs_vec += implicit_rhs;

   // Reduce system
   auto K_hat_red = Timer::measure(
       [&]{ return mattools::reduce_mat_to_mat(K_hat, full_lhs_vec); },
       "reducing K_hat to solver-ready form" );

   auto rhs_red = Timer::measure(
       [&]{ return mattools::reduce_vec_to_vec(full_rhs_vec, full_lhs_vec); },
       "reducing load vector to solver-ready form" );

   K_hat_red.makeCompressed();

   // -------------------------------------------------------------------------
   // 5) Solve static preload displacements
   // -------------------------------------------------------------------------
   auto sol = Timer::measure(
       [&]{ return solve(device, method, K_hat_red, rhs_red); },
       "solving linear system for pre-stress" );

   DynamicVector u_active = sol.head(sol.rows() - n);
   DynamicVector u_full   = mattools::expand_vec_to_vec(u_active, active_lhs_vec);
   DynamicMatrix U_mat    = mattools::expand_vec_to_mat(active_dof_idx_mat, u_full);

   // -------------------------------------------------------------------------
   // 6) Compute integration point stresses
   // -------------------------------------------------------------------------
   IPData ip_stress;
   {
       NodeData U_row = U_mat;
       IPData ip_strain_unused;
       std::tie(ip_stress, ip_strain_unused) =
           Timer::measure([&]{ return m_model->compute_ip_stress_strain(U_row); },
                          "computing IP stress/strain for Kg");
   }

   // -------------------------------------------------------------------------
   // 7) Assemble geometric stiffness Kg
   // -------------------------------------------------------------------------
   auto Kg_act = Timer::measure(
       [&]{ return m_model->build_geom_stiffness_matrix(active_dof_idx_mat, ip_stress); },
       "assembling active geometric stiffness matrix" );

   auto Kg_hat = Timer::measure(
       [&]() {
           SparseMatrix G(m + n, m + n);
           TripletList T; T.reserve(Kg_act.nonZeros());
           for (int k = 0; k < Kg_act.outerSize(); ++k)
               for (SparseMatrix::InnerIterator it(Kg_act, k); it; ++it)
                   T.emplace_back(it.row(), it.col(), it.value());
           G.setFromTriplets(T.begin(), T.end());
           return G;
       },
       "assembling Kghat" );

   write_mtx<Precision>("K_norm.mtx", K_hat , 0.0001);
   write_mtx<Precision>("K_geom.mtx", Kg_hat, 0.0001);

   auto Kg_hat_red = Timer::measure(
       [&]{ return mattools::reduce_mat_to_mat(Kg_hat, full_lhs_vec); },
       "reducing Kghat to solver-ready form" );

   K_hat_red.makeCompressed();
   Kg_hat_red.makeCompressed();

   // -------------------------------------------------------------------------
   // 8) Estimate first eigenvalue and choose shift σ
   // -------------------------------------------------------------------------
   double lambda_est = Timer::measure(
       [&]{ return estimate_first_eigenvalue(K_hat_red, Kg_hat_red, u_active); },
       "estimating first eigenvalue (Rayleigh quotient)" );

   double sigma_val = lambda_est / 1e6;   // scale down strongly
    // if (sigma_val <= 0.0) sigma_val *= -1;
   if (!std::isfinite(sigma_val))
       sigma_val = 1e-9; // fallback safe value

   // -------------------------------------------------------------------------
   // 9) Solve buckling eigenvalue problem
   // -------------------------------------------------------------------------
   SparseMatrix B = -Kg_hat_red; // ensure positive buckling factors

   const Index active_system_dofs = rhs_red.rows() - n;
   const Index neigs = std::min<Index>(num_eigenvalues,
                                       std::max<Index>(1, active_system_dofs));

   solver::EigenOpts opts;
   opts.mode  = solver::EigenMode::Buckling;
   opts.sigma = sigma_val;
   opts.sort  = solver::EigenOpts::Sort::LargestMagn;

   const int ncv = std::min<int>(static_cast<int>(K_hat_red.rows()),
                                 std::max<int>(3 * static_cast<int>(neigs) + 20,
                                               static_cast<int>(neigs) + 2));
   opts.ncv = ncv;

   auto eig_pairs = Timer::measure(
       [&]{ return solver::eigs(solver::CPU, K_hat_red, B,
                                 static_cast<int>(neigs), opts); },
       "solving generalized EVP for buckling"
   );

   auto modes = make_and_sort_modes(eig_pairs);
   for (auto& m : modes)
       m.expand_mode_shape(n, active_lhs_vec, active_dof_idx_mat);

   display_buckling_results(modes, opts.sigma, ncv, static_cast<int>(neigs));

   // -------------------------------------------------------------------------
   // 10) Write results
   // -------------------------------------------------------------------------
   m_writer->add_loadcase(m_id);
   {
       DynamicVector lambdas(modes.size());
       for (size_t i = 0; i < modes.size(); ++i) {
           lambdas(i) = modes[i].lambda;
           m_writer->write_eigen_matrix(modes[i].mode_shape_mat,
                                        "BUCKLING_MODE_" + std::to_string(i+1));
       }
       m_writer->write_eigen_matrix(DynamicMatrix(lambdas), "BUCKLING_FACTORS");
   }

   logging::info(true, "Buckling analysis completed.");
}

}} // namespace fem::loadcase
