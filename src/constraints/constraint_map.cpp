/******************************************************************************
 * @file constraint_map.cpp
 * @brief Implementation of the immutable null-space map `u = u_p + T q`.
 *
 * -----------------------------------------------------------------------------
 * ## What this file provides
 *
 * Efficient, allocation-light implementations of:
 *   • Applying the affine map and its transpose:
 *         u = u_p + T q,         z = Tᵀ y
 *   • Assembling reduced operators and RHS:
 *         A = Tᵀ K T,            B = Tᵀ K_g T,         b = Tᵀ (f − K u_p)
 *   • Matrix-free operator wrappers OpA/OpB for iterative solvers.
 *   • Recovery of the full solution and computation of reaction forces:
 *         u = u_p + T q,         r = K u − f
 *
 * The map is constructed by ConstraintBuilder and then treated as immutable.
 *
 * -----------------------------------------------------------------------------
 * ## Why this matters (ELI5)
 *
 * Constraints `C u = d` are “strings” tying some DOFs (slaves) to others
 * (masters). After building this map, you only solve for the free knobs `q`
 * (masters); the rest are computed automatically via `T` and `u_p`. This
 * works for:
 *   • Linear static:      Tᵀ K T q = Tᵀ (f − K u_p)
 *   • Modal (frequency):  Tᵀ K T q = λ Tᵀ M T q
 *   • Buckling:           Tᵀ K T q = λ Tᵀ K_g T q
 *
 * -----------------------------------------------------------------------------
 * ## Performance notes
 *
 * - `apply_T` / `apply_Tt` are O(nm + nnz(X)), with nm = n_master().
 * - OpA/OpB perform one (sparse) K·(T·x) matvec and a cheap Tᵀ application.
 * - Assembling A/B forms explicit reduced matrices; for very large n, prefer
 *   matrix-free via OpA/OpB and iterative methods (e.g., CG).
 *
 * @see constraint_map.h for the class interface and usage example.
 ******************************************************************************/

#include "constraint_map.h"
#include "../core/logging.h"

namespace fem::constraint {

/******************************************************************************
 * @brief Apply full map: u = u_p + T q.
 *
 * The sparse `T` is not explicitly multiplied; we exploit the structure:
 *   - For master rows, `u[masters[k]] += q[k]`.
 *   - For slave rows, we accumulate `(X row i) · q` into `u[slaves[i]]`.
 *
 * @param q Reduced vector (size = n_master()).
 * @param u Output full vector (size = n_full()).
 ******************************************************************************/
void ConstraintMap::apply_T(const DynamicVector& q, DynamicVector& u) const {
    // Basic size sanity (runtime defensive checks without hard aborts).
    if (q.size() != nm_) {
        logging::warning(false, "[ConstraintMap::apply_T] q has size ", q.size(),
                         " but nm_ = ", nm_, ". Result may be invalid.");
    }

    u = u_p_; // start from particular solution (zeros for homogeneous constraints)

    // Masters: identity injection
    for (int k = 0; k < static_cast<int>(masters_.size()); ++k) {
        u[masters_[k]] += q[k];
    }

    // Slaves: add X * q into the slave rows
    for (int j = 0; j < X_.outerSize(); ++j) {
        for (SparseMatrix::InnerIterator it(X_, j); it; ++it) {
            u[slaves_[it.row()]] += it.value() * q[j];
        }
    }
}

/******************************************************************************
 * @brief Apply transpose map: z = Tᵀ y.
 *
 * Using the same structure without forming dense matrices:
 *   - Master part: z[k] accumulates y at the master DOF row.
 *   - Slave part:  z[j] accumulates (Xᵀ y_slaves)[j].
 *
 * @param y Full vector (size = n_full()).
 * @param z Output reduced vector (size = n_master()).
 ******************************************************************************/
void ConstraintMap::apply_Tt(const DynamicVector& y, DynamicVector& z) const {
    if (y.size() != n_) {
        logging::warning(false, "[ConstraintMap::apply_Tt] y has size ", y.size(),
                         " but n_ = ", n_, ". Result may be invalid.");
    }

    z.setZero(nm_);

    // Masters: identity gather
    for (int k = 0; k < static_cast<int>(masters_.size()); ++k) {
        z[k] += y[masters_[k]];
    }

    // Slaves: add Xᵀ * y_slaves
    for (int j = 0; j < X_.outerSize(); ++j) {
        Precision acc = 0;
        for (SparseMatrix::InnerIterator it(X_, j); it; ++it) {
            acc += it.value() * y[slaves_[it.row()]];
        }
        z[j] += acc;
    }
}

/******************************************************************************
 * @brief Convenience overload: return u = u_p + T q by value.
 ******************************************************************************/
DynamicVector ConstraintMap::apply_T(const DynamicVector& q) const {
    DynamicVector u(n_);
    apply_T(q, u);
    return u;
}

/******************************************************************************
 * @brief Convenience overload: return z = Tᵀ y by value.
 ******************************************************************************/
DynamicVector ConstraintMap::apply_Tt(const DynamicVector& y) const {
    DynamicVector z(nm_);
    apply_Tt(y, z);
    return z;
}

/******************************************************************************
 * @brief Assemble reduced operator A = Tᵀ K T explicitly.
 *
 * Prefer the matrix-free OpA for very large problems. This version performs
 * two sparse multiplications: (K * T) and (Tᵀ * (K * T)).
 *
 * @param K Full (typically SPD) matrix.
 * @return Reduced matrix A (size nm × nm).
 ******************************************************************************/
SparseMatrix ConstraintMap::assemble_A(const SparseMatrix& K) const {
    SparseMatrix KT = K * T_;
    SparseMatrix A  = T_.transpose() * KT;
    A.makeCompressed();
    return A;
}

/******************************************************************************
 * @brief Assemble reduced operator B = Tᵀ K_g T explicitly.
 *
 * Useful for generalized EVP problems (e.g., frequency or buckling).
 *
 * @param Kg Full geometric or mass matrix.
 * @return Reduced matrix B (size nm × nm).
 ******************************************************************************/
SparseMatrix ConstraintMap::assemble_B(const SparseMatrix& Kg) const {
    SparseMatrix KgT = Kg * T_;
    SparseMatrix B   = T_.transpose() * KgT;
    B.makeCompressed();
    return B;
}

/******************************************************************************
 * @brief Assemble reduced RHS b = Tᵀ (f − K u_p).
 *
 * @param K Full matrix K.
 * @param f Full load vector f.
 * @return Reduced RHS b (size nm).
 ******************************************************************************/
DynamicVector ConstraintMap::assemble_b(const SparseMatrix& K, const DynamicVector& f) const {
    if (f.size() != n_) {
        logging::warning(false, "[ConstraintMap::assemble_b] f has size ", f.size(),
                         " but n_ = ", n_, ". Result may be invalid.");
    }
    DynamicVector tmp = f - K * u_p_;
    return apply_Tt(tmp);
}

/******************************************************************************
 * @class ConstraintMap::OpA
 * @brief Matrix-free operator: y = Tᵀ (K (T x)).
 *
 * Reuses internal buffers to minimize allocations across repeated calls.
 ******************************************************************************/
ConstraintMap::OpA::OpA(const ConstraintMap& map, const SparseMatrix& K)
    : map_(map),
      K_(K),
      u_full_(DynamicVector::Zero(map.n_)),
      y_full_(DynamicVector::Zero(map.n_)) {}

/******************************************************************************
 * @brief Perform y = Tᵀ K T x.
 ******************************************************************************/
void ConstraintMap::OpA::perform_op(const DynamicVector& x, DynamicVector& y) const {
    map_.apply_T(x, u_full_);
    y_full_ = K_ * u_full_;
    map_.apply_Tt(y_full_, y);
}

/******************************************************************************
 * @class ConstraintMap::OpB
 * @brief Matrix-free operator: y = Tᵀ (K_g (T x)).
 ******************************************************************************/
ConstraintMap::OpB::OpB(const ConstraintMap& map, const SparseMatrix& Kg)
    : map_(map),
      Kg_(Kg),
      u_full_(DynamicVector::Zero(map.n_)),
      y_full_(DynamicVector::Zero(map.n_)) {}

/******************************************************************************
 * @brief Perform y = Tᵀ K_g T x.
 ******************************************************************************/
void ConstraintMap::OpB::perform_op(const DynamicVector& x, DynamicVector& y) const {
    map_.apply_T(x, u_full_);
    y_full_ = Kg_ * u_full_;
    map_.apply_Tt(y_full_, y);
}

/******************************************************************************
 * @brief Recover full vector u = u_p + T q.
 ******************************************************************************/
DynamicVector ConstraintMap::recover_u(const DynamicVector& q) const {
    return apply_T(q);
}

/******************************************************************************
 * @brief Compute full residual / reactions r = K u − f.
 *
 * For a correctly reduced (constrained) solve, the projected residual must be
 * (near) zero:  ‖Tᵀ r‖ ≈ 0. The residual then lies in range(Cᵀ) and represents
 * **reaction forces/moments** at constrained DOFs.
 *
 * @param K Full matrix K.
 * @param f Full RHS f.
 * @param q Reduced solution q.
 * @return r = K u − f (size n).
 ******************************************************************************/
DynamicVector ConstraintMap::reactions(const SparseMatrix& K,
                                       const DynamicVector& f,
                                       const DynamicVector& q) const {
    DynamicVector u = recover_u(q);
    return K * u - f;
}

} // namespace fem::constraint
