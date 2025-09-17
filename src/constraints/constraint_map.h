/******************************************************************************
 * @file constraint_map.h
 * @brief Immutable null-space map `u = u_p + T q` with fast application helpers.
 *
 * @date    14.09.2025 (API stable; sparse-only apply/assemble paths)  @author Finn
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include <vector>

namespace fem::constraint {

class ConstraintMap {
public:
    ConstraintMap() = default;

    // sizes
    Index n_full()   const { return n_;  }
    Index n_master() const { return nm_; }
    Index n_slave()  const { return r_;  }

    // indices
    const std::vector<Index>& master_idx() const { return masters_; }
    const std::vector<Index>& slave_idx()  const { return slaves_;  }

    // components
    const DynamicVector& u_p() const { return u_p_; }
    const SparseMatrix&  X()   const { return X_;   }
    const SparseMatrix&  T()   const { return T_;   }

    // application
    void          apply_T (const DynamicVector& q, DynamicVector& u) const;
    void          apply_Tt(const DynamicVector& y, DynamicVector& z) const;
    DynamicVector apply_T (const DynamicVector& q) const;
    DynamicVector apply_Tt(const DynamicVector& y) const;

    // assembled reduced operators
    SparseMatrix  assemble_A(const SparseMatrix& K)  const; // A = Tᵀ K T
    SparseMatrix  assemble_B(const SparseMatrix& Kg) const; // B = Tᵀ Kg T
    DynamicVector assemble_b(const SparseMatrix& K, const DynamicVector& f) const;

    // matrix-free ops
    class OpA {
    public:
        OpA(const ConstraintMap& map, const SparseMatrix& K);
        void perform_op(const DynamicVector& x, DynamicVector& y) const; // y = Tᵀ K T x
    private:
        const ConstraintMap& map_;
        const SparseMatrix&  K_;
        mutable DynamicVector u_full_, y_full_;
    };
    class OpB {
    public:
        OpB(const ConstraintMap& map, const SparseMatrix& Kg);
        void perform_op(const DynamicVector& x, DynamicVector& y) const; // y = Tᵀ Kg T x
    private:
        const ConstraintMap& map_;
        const SparseMatrix&  Kg_;
        mutable DynamicVector u_full_, y_full_;
    };

    // recovery / reactions
    DynamicVector recover_u(const DynamicVector& q) const;
    DynamicVector reactions(const SparseMatrix& K, const DynamicVector& f, const DynamicVector& q) const;

public: // constructed by ConstraintBuilder
    Index n_{0};
    Index r_{0};
    Index nm_{0};

    std::vector<Index> masters_;
    std::vector<Index> slaves_;

    SparseMatrix  X_;
    DynamicVector u_p_;
    SparseMatrix  T_;
};

} // namespace fem::constraint
