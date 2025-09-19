/******************************************************************************
 * @file constraint_map.cpp
 * @brief Implementation of the immutable null-space map `u = u_p + T q`.
 ******************************************************************************/

#include "constraint_map.h"

#include "../core/logging.h"

#include <cmath>

namespace fem::constraint {

void ConstraintMap::apply_T(const DynamicVector& q, DynamicVector& u) const {
    if (q.size() != nm_) {
        logging::warning(false, "[ConstraintMap::apply_T] q.size()=", q.size(), " nm_=", nm_);
    }
    u = u_p_; // start with particular solution (often zero)

    // Masters: identity add
    for (int k = 0; k < (int)masters_.size(); ++k) {
        u[masters_[k]] += q[k];
    }
    // Slaves: add X * q into slave rows
    for (int j = 0; j < X_.outerSize(); ++j) {
        for (SparseMatrix::InnerIterator it(X_, j); it; ++it) {
            u[slaves_[it.row()]] += it.value() * q[j];
        }
    }
}

void ConstraintMap::apply_Tt(const DynamicVector& y, DynamicVector& z) const {
    if (y.size() != n_) {
        logging::warning(false, "[ConstraintMap::apply_Tt] y.size()=", y.size(), " n_=", n_);
    }
    z.setZero(nm_);

    // Masters: gather identity
    for (int k = 0; k < (int)masters_.size(); ++k) {
        z[k] += y[masters_[k]];
    }
    // Slaves: z += Xᵀ * y_slaves
    for (int j = 0; j < X_.outerSize(); ++j) {
        Precision acc = 0;
        for (SparseMatrix::InnerIterator it(X_, j); it; ++it) {
            acc += it.value() * y[slaves_[it.row()]];
        }
        z[j] += acc;
    }
}

DynamicVector ConstraintMap::apply_T(const DynamicVector& q) const {
    DynamicVector u(n_);
    apply_T(q, u);
    return u;
}
DynamicVector ConstraintMap::apply_Tt(const DynamicVector& y) const {
    DynamicVector z(nm_);
    apply_Tt(y, z);
    return z;
}

SparseMatrix ConstraintMap::assemble_A(const SparseMatrix& K) const {
    // A = Tᵀ (K T) ; avoid dense intermediates
    SparseMatrix KT = K * T_;
    SparseMatrix A  = T_.transpose() * KT;
    A.makeCompressed();
    return A;
}

SparseMatrix ConstraintMap::assemble_B(const SparseMatrix& Kg) const {
    SparseMatrix KgT = Kg * T_;
    SparseMatrix B   = T_.transpose() * KgT;
    B.makeCompressed();
    return B;
}

DynamicVector ConstraintMap::assemble_b(const SparseMatrix& K, const DynamicVector& f) const {
    DynamicVector Ku = K * u_p_;
    DynamicVector tmp = f - Ku;
    return apply_Tt(tmp);
}

ConstraintMap::OpA::OpA(const ConstraintMap& map, const SparseMatrix& K)
    : map_(map), K_(K),
      u_full_(DynamicVector::Zero(map.n_)),
      y_full_(DynamicVector::Zero(map.n_)) {}

void ConstraintMap::OpA::perform_op(const DynamicVector& x, DynamicVector& y) const {
    map_.apply_T(x, u_full_);      // u = T x
    y_full_ = K_ * u_full_;        // y_full = K u
    map_.apply_Tt(y_full_, y);     // y = Tᵀ y_full
}

ConstraintMap::OpB::OpB(const ConstraintMap& map, const SparseMatrix& Kg)
    : map_(map), Kg_(Kg),
      u_full_(DynamicVector::Zero(map.n_)),
      y_full_(DynamicVector::Zero(map.n_)) {}

void ConstraintMap::OpB::perform_op(const DynamicVector& x, DynamicVector& y) const {
    map_.apply_T(x, u_full_);
    y_full_ = Kg_ * u_full_;
    map_.apply_Tt(y_full_, y);
}

DynamicVector ConstraintMap::recover_u(const DynamicVector& q) const {
    return apply_T(q);
}

DynamicVector ConstraintMap::reactions(const SparseMatrix& K,
                                       const DynamicVector& f,
                                       const DynamicVector& q) const {
    DynamicVector u = recover_u(q);
    return K * u - f;
}

} // namespace fem::constraint
