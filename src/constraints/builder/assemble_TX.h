/******************************************************************************
* @file assemble_TX.h
 * @brief Assemble T (sparse) and X (sparse) from X_cols and partition info.
 ******************************************************************************/
#pragma once
#include "../../core/types_eig.h"
#include "build_x.h"
#include <vector>

namespace fem { namespace constraint {

struct AssembleInput {
    int n;                              ///< total DOFs
    int r;                              ///< rank
    std::vector<int>      slaves_loc;   ///< local → used
    std::vector<int>      masters_loc;  ///< local → used
    std::vector<int>      used;         ///< used → global
    std::vector<Index>    masters_glob; ///< final master order
};

struct AssembleOutput {
    SparseMatrix T;   ///< n x nm
    SparseMatrix X;   ///< r x nm (zeros on trailing never-used masters)
};

AssembleOutput assemble_T_and_X(const AssembleInput& in, const XCols& Xc);

}} // namespace
