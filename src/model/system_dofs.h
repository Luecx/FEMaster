/**
 * @file system_dofs.h
 * @brief Declares the node-by-component system DOF numbering.
 *
 * `SystemDofs` owns the activation mask and the contiguous global numbering of
 * the algebraic system. It is built in two phases: components are activated
 * first, then `numerate()` assigns zero-based indices in node-major order.
 * After numbering, it maps nodal fields to system vectors and back and builds
 * the sparsity pattern of the global matrices.
 *
 * @see src/model/system_dofs.cpp
 * @author Finn Eggers
 * @date 19.11.2026
 */

#pragma once

#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../data/field.h"

#include <vector>

namespace fem::model {

class SystemDofs {
public:
    // Empty numbering without nodes; used where no system mapping is required
    SystemDofs() = default;

    // All components of all nodes start inactive
    SystemDofs(Index node_count, Index dof_count);

    // Activation, only valid before numerate()
    void activate(Index node, Index dof);
    void activate(Index node, const Dofs& dofs);

    // Assigns contiguous indices to active components in node-major order
    void numerate();

    // Global index of a component, or -1 if inactive. Unchecked for assembly loops.
    int operator()(Index node, Index dof) const { return _ids(node, dof); }

    bool  active   (Index node, Index dof) const { return _mask(node, dof); }
    bool  numerated()                      const { return _numerated; }
    Index size     ()                      const { return _size; }
    Index nodes    ()                      const { return _mask.rows(); }
    Index width    ()                      const { return _mask.cols(); }

    // Maps between nodal fields of size nodes() x width() and system vectors of size size()
    DynamicVector reduce(const Field&         field ) const;
    Field         expand(const DynamicVector& vector) const;

    // Zero-valued size() x size() matrix with every coupling admitted by the element topology
    SparseMatrix  pattern(const std::vector<ElementPtr>& elements) const;

private:
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _mask;
    Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _ids;
    Index                                                                _size      = 0;
    bool                                                                 _numerated = false;
};

} // namespace fem::model
