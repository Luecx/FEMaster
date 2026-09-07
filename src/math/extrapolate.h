/**
 * @file extrapolate.h
 * @brief Declares polynomial extrapolation between points in natural coordinates.
 *
 * The module builds a linear recovery operator from values known at source
 * points to arbitrary target points. A caller explicitly selects the polynomial
 * basis used for the reconstruction; element topology and quadrature remain
 * responsibilities of the calling formulation.
 *
 * @author Finn Eggers
 * @date 07.09.2026
 */

#pragma once

#include "../core/types_eig.h"

#include <initializer_list>

namespace fem::math {

/**
 * @brief Polynomial basis functions available for reference-space recovery.
 *
 * The symbols use the natural coordinates `(r, s, t)`. Higher-order entries are
 * simple monomials and may be combined freely by the caller.
 */
enum class ExtrapolationBasis {
    F1,

    FR,
    FS,
    FT,

    FRR,
    FSS,
    FTT,

    FRS,
    FRT,
    FST,
    FRST,

    FRRS,
    FRRT,
    FSSR,
    FSST,
    FTTR,
    FTTS,

    FRRST,
    FRSST,
    FRSTT,
};

RowMatrix extrapolate(const RowMatrix&                              source_points,
                      const RowMatrix&                              target_points,
                      std::initializer_list<ExtrapolationBasis>     basis);

} // namespace fem::math
