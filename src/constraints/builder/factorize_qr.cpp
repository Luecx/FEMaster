/******************************************************************************
 * @file factorize_qr.cpp
 * @brief Wraps Eigen sparse QR factorisation with custom settings.
 *
 * The implementation forwards the configured pivot threshold and exposes the
 * resulting `R` factor and column permutation.
 *
 * @see src/constraints/builder/factorize_qr.h
 * @see src/constraints/builder/preprocess.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "factorize_qr.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

namespace fem {
namespace constraint {

/******************************************************************************
 * @copydoc factorize_sparse_qr
 ******************************************************************************/
void factorize_sparse_qr(const SparseMatrix& C_use, const QRSettings& settings, QRResult& out) {
    out.qr.setPivotThreshold(settings.pivot_rel <= 0 ? 0 : std::min<Precision>(0.1, settings.pivot_rel));
    out.qr.compute(C_use);
    logging::error(out.qr.info() == Eigen::Success, "[ConstraintBuilder] SparseQR failed.");

    out.R = out.qr.matrixR();
    out.P = out.qr.colsPermutation();
}

} // namespace constraint
} // namespace fem
