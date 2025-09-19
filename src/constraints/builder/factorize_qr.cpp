/******************************************************************************
* @file factorize_qr.cpp
* @brief SparseQR factorization wrapper with configurable pivot threshold.
******************************************************************************/

#include "factorize_qr.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

namespace fem { namespace constraint {

void factorize_sparse_qr(const SparseMatrix& C_use, const QRSettings& s, QRResult& out)
{
    out.qr.setPivotThreshold(s.pivot_rel <= 0 ? 0 : std::min<Precision>(0.1, s.pivot_rel));
    out.qr.compute(C_use);
    logging::error(out.qr.info() == Eigen::Success, "[ConstraintBuilder] SparseQR failed.");

    out.R = out.qr.matrixR();
    out.P = out.qr.colsPermutation();
}

}} // namespace fem::constraint
