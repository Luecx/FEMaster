/******************************************************************************
 * @file constraint_map.cpp
 * @brief Implements the null-space map utilities used by constraint handling.
 *
 * The functions defined here operate on the transformation `u = u_p + T q` that
 * eliminates constrained DOFs from the global system.
 *
 * @see src/constraints/constraint_map.h
 * @see src/constraints/builder/builder.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "constraint_map.h"

#include "../core/logging.h"

#include <cmath>

namespace fem {
namespace constraint {

/******************************************************************************
 * @copydoc ConstraintMap::apply_T(const DynamicVector&,DynamicVector&) const
 ******************************************************************************/
void ConstraintMap::apply_T(const DynamicVector& q, DynamicVector& u) const {
    if (q.size() != nm_) {
        logging::warning(false, "[ConstraintMap::apply_T] q.size()=", q.size(), " nm_=", nm_);
    }
    u = u_p_;

    for (int k = 0; k < static_cast<int>(masters_.size()); ++k) {
        u[masters_[k]] += q[k];
    }
    for (int j = 0; j < X_.outerSize(); ++j) {
        for (SparseMatrix::InnerIterator it(X_, j); it; ++it) {
            u[slaves_[it.row()]] += it.value() * q[j];
        }
    }
}

/******************************************************************************
 * @copydoc ConstraintMap::apply_Tt(const DynamicVector&,DynamicVector&) const
 ******************************************************************************/
void ConstraintMap::apply_Tt(const DynamicVector& y, DynamicVector& z) const {
    if (y.size() != n_) {
        logging::warning(false, "[ConstraintMap::apply_Tt] y.size()=", y.size(), " n_=", n_);
    }
    z.setZero(nm_);

    for (int k = 0; k < static_cast<int>(masters_.size()); ++k) {
        z[k] += y[masters_[k]];
    }
    for (int j = 0; j < X_.outerSize(); ++j) {
        Precision acc = 0;
        for (SparseMatrix::InnerIterator it(X_, j); it; ++it) {
            acc += it.value() * y[slaves_[it.row()]];
        }
        z[j] += acc;
    }
}

/******************************************************************************
 * @copydoc ConstraintMap::apply_T(const DynamicVector&) const
 ******************************************************************************/
DynamicVector ConstraintMap::apply_T(const DynamicVector& q) const {
    DynamicVector u(n_);
    apply_T(q, u);
    return u;
}

/******************************************************************************
 * @copydoc ConstraintMap::apply_Tt(const DynamicVector&) const
 ******************************************************************************/
DynamicVector ConstraintMap::apply_Tt(const DynamicVector& y) const {
    DynamicVector z(nm_);
    apply_Tt(y, z);
    return z;
}

/******************************************************************************
 * @copydoc ConstraintMap::assemble_A(const SparseMatrix&) const
 ******************************************************************************/
SparseMatrix ConstraintMap::assemble_A(const SparseMatrix& K) const {
    SparseMatrix KT = K * T_;
    SparseMatrix A = T_.transpose() * KT;
    A.makeCompressed();
    return A;
}

/******************************************************************************
 * @copydoc ConstraintMap::assemble_B(const SparseMatrix&) const
 ******************************************************************************/
SparseMatrix ConstraintMap::assemble_B(const SparseMatrix& Kg) const {
    SparseMatrix KgT = Kg * T_;
    SparseMatrix B = T_.transpose() * KgT;
    B.makeCompressed();
    return B;
}

/******************************************************************************
 * @copydoc ConstraintMap::assemble_b(const SparseMatrix&,const DynamicVector&) const
 ******************************************************************************/
DynamicVector ConstraintMap::assemble_b(const SparseMatrix& K, const DynamicVector& f) const {
    DynamicVector Ku = K * u_p_;
    DynamicVector tmp = f - Ku;
    return apply_Tt(tmp);
}

/******************************************************************************
 * @copydoc ConstraintMap::OpA::OpA
 ******************************************************************************/
ConstraintMap::OpA::OpA(const ConstraintMap& map, const SparseMatrix& K)
    : map_(map), K_(K),
      u_full_(DynamicVector::Zero(map.n_)),
      y_full_(DynamicVector::Zero(map.n_)) {}

/******************************************************************************
 * @copydoc ConstraintMap::OpA::perform_op
 ******************************************************************************/
void ConstraintMap::OpA::perform_op(const DynamicVector& x, DynamicVector& y) const {
    map_.apply_T(x, u_full_);
    y_full_ = K_ * u_full_;
    map_.apply_Tt(y_full_, y);
}

/******************************************************************************
 * @copydoc ConstraintMap::OpB::OpB
 ******************************************************************************/
ConstraintMap::OpB::OpB(const ConstraintMap& map, const SparseMatrix& Kg)
    : map_(map), Kg_(Kg),
      u_full_(DynamicVector::Zero(map.n_)),
      y_full_(DynamicVector::Zero(map.n_)) {}

/******************************************************************************
 * @copydoc ConstraintMap::OpB::perform_op
 ******************************************************************************/
void ConstraintMap::OpB::perform_op(const DynamicVector& x, DynamicVector& y) const {
    map_.apply_T(x, u_full_);
    y_full_ = Kg_ * u_full_;
    map_.apply_Tt(y_full_, y);
}

/******************************************************************************
 * @copydoc ConstraintMap::recover_u
 ******************************************************************************/
DynamicVector ConstraintMap::recover_u(const DynamicVector& q) const {
    return apply_T(q);
}

/******************************************************************************
 * @copydoc ConstraintMap::reactions
 ******************************************************************************/
DynamicVector ConstraintMap::reactions(const SparseMatrix& K,
                                       const DynamicVector& f,
                                       const DynamicVector& q) const {
    DynamicVector u = recover_u(q);
    return K * u - f;
}

} // namespace constraint
} // namespace fem
