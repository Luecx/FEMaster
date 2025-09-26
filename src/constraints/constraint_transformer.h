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

#include <iostream>

namespace fem::constraint {

class ConstraintTransformer {
public:
    struct BuildOptions {
        ConstraintSet::Options     set;
        ConstraintBuilder::Options builder;
    };

    ConstraintTransformer(const Equations& eqs,
                          const SystemDofIds& dofs,
                          Index n_dofs,
                          const BuildOptions& opt = {})
    {
        set_.equations = eqs;
        set_.opt = opt.set;
        std::cout << "assembling set..." << std::endl;
        set_.assemble(dofs, n_dofs);
        std::cout << "building map..." << std::endl;
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

private:
    ConstraintSet             set_;
    ConstraintMap             map_;
    ConstraintBuilder::Report report_;
};

} // namespace fem::constraint
