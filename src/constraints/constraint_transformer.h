/******************************************************************************
 * @file constraint_transformer.h
 * @brief One-stop façade to build and use the constraint null-space map.
 *
 * @date    14.09.2025                                                    @author Finn
 ******************************************************************************/

#pragma once

#include "builder/builder.h"
#include "constraint_map.h"
#include "constraint_set.h"
#include <limits>
#include <iostream>

namespace fem::constraint {



class ConstraintTransformer {
public:
    struct BuildOptions {
        ConstraintSet::Options     set;
        ConstraintBuilder::Options builder;
    };


    struct StaticCheckOptions {
        // Relative tolerances
        Precision tol_constraint_rel = 1e-10; // ||C u - d|| / max(1, ||d||)
        Precision tol_reduced_rel    = 1e-8;  // ||T^T (K u - f)|| / denom_red
        // Optional: full residual tolerance (set inf to ignore)
        Precision tol_full_rel       = std::numeric_limits<Precision>::infinity();

        // Denominator for reduced relative residual:
        // denom_red = max(1, ||T^T f||, ||T^T K u||) for robustness
    };

    struct StaticCheckResult {
        // Raw norms
        Precision norm_u        = 0;
        Precision norm_Cu_d     = 0;
        Precision norm_resid    = 0; // ||K u - f||
        Precision norm_reduced  = 0; // ||T^T (K u - f)||

        // Relatives
        Precision rel_Cu_d      = 0;
        Precision rel_resid     = 0; // full-space rel residual
        Precision rel_reduced   = 0; // reduced-space rel residual

        // Denominators used (handy for logs)
        Precision denom_Cu_d    = 1;
        Precision denom_full    = 1;
        Precision denom_reduced = 1;

        // Pass flags
        bool pass_constraints   = false;
        bool pass_reduced_eq    = false;
        bool pass_full_resid    = true; // ignored if tol_full_rel is inf
    };

    ConstraintTransformer(const Equations& eqs,
                          const SystemDofIds& dofs,
                          Index n_dofs,
                          const BuildOptions& opt = {})
    {
        set_.equations = eqs;
        set_.opt = opt.set;
        set_.assemble(dofs, n_dofs);
        std::tie(map_, report_) = ConstraintBuilder::build(set_, opt.builder);
    }

    // Access
    const ConstraintMap&             map()    const { return map_;    }
    const ConstraintSet&             set()    const { return set_;    }
    const ConstraintBuilder::Report& report() const { return report_; }

    // Quick
    bool  homogeneous() const { return report_.homogeneous; }
    bool  feasible()    const { return report_.feasible; }
    Index rank()        const { return report_.rank; }
    Index n_master()    const { return map_.n_master(); }

    // Reduced assemblies
    SparseMatrix  assemble_A(const SparseMatrix& K)  const { return map_.assemble_A(K);  }
    SparseMatrix  assemble_B(const SparseMatrix& Kg) const { return map_.assemble_B(Kg); }
    DynamicVector assemble_b(const SparseMatrix& K, const DynamicVector& f) const {
        return map_.assemble_b(K, f);
    }

    // Matrix-free operators
    ConstraintMap::OpA opA(const SparseMatrix& K)  const { return ConstraintMap::OpA(map_, K);  }
    ConstraintMap::OpB opB(const SparseMatrix& Kg) const { return ConstraintMap::OpB(map_, Kg); }

    // Applications & utilities
    void          apply_T (const DynamicVector& q, DynamicVector& u) const { map_.apply_T(q, u); }
    void          apply_Tt(const DynamicVector& y, DynamicVector& z) const { map_.apply_Tt(y, z); }
    DynamicVector recover_u(const DynamicVector& q) const { return map_.recover_u(q); }
    DynamicVector reactions(const SparseMatrix& K, const DynamicVector& f, const DynamicVector& q) const {
        return map_.reactions(K, f, q);
    }


    // Main implementations (no defaults)
    StaticCheckResult check_static(const SparseMatrix& K,
                                   const DynamicVector& f,
                                   const DynamicVector& u,
                                   StaticCheckOptions opt) const;

    void print_checklist(const StaticCheckResult& r,
                         StaticCheckOptions opt) const;

    // Convenience overloads (use default options)
    StaticCheckResult check_static(const SparseMatrix& K,
                                   const DynamicVector& f,
                                   const DynamicVector& u) const;

    void print_checklist(const StaticCheckResult& r) const;

private:
    ConstraintSet             set_;
    ConstraintMap             map_;
    ConstraintBuilder::Report report_;
};

} // namespace fem::constraint
