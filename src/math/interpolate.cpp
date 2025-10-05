/**
 * @file interpolate.cpp
 * @brief Implements polynomial interpolation for scattered nodal data.
 *
 * @see src/math/interpolate.h
 */

#include "interpolate.h"

#include <Eigen/Eigen>
#include <Eigen/QR>

#include <cmath>

namespace fem {
namespace math {
namespace interpolate {

namespace {

/**
 * @brief Returns the number of coefficients required for interpolation order `F`.
 */
template<InterpolationFunction F>
constexpr int term_count() {
    if constexpr (F == CONSTANT) {
        return 1;
    } else if constexpr (F == LINEAR) {
        return 4;
    } else if constexpr (F == BILINEAR) {
        return 7;
    } else if constexpr (F == QUADRATIC) {
        return 10;
    } else if constexpr (F == BILINQUAD) {
        return 17;
    } else if constexpr (F == CUBIC) {
        return 20;
    } else {
        return 0;
    }
}

/**
 * @brief Demotes the interpolation order to the next simpler option.
 */
template<InterpolationFunction F>
constexpr InterpolationFunction demote_order() {
    if constexpr (F == CUBIC) {
        return BILINQUAD;
    } else if constexpr (F == BILINQUAD) {
        return QUADRATIC;
    } else if constexpr (F == QUADRATIC) {
        return BILINEAR;
    } else if constexpr (F == BILINEAR) {
        return LINEAR;
    } else if constexpr (F == LINEAR) {
        return CONSTANT;
    } else {
        return CONSTANT;
    }
}

/**
 * @brief Evaluates the polynomial defined by `coefficients` at `position`.
 */
template<InterpolationFunction F>
Precision evaluate(const DynamicVector& coefficients, const Vec3& position) {
    Precision result = 0.0;
    int idx = 0;

    if constexpr (F >= CONSTANT) {
        result += coefficients(idx++);
    }
    if constexpr (F >= LINEAR) {
        result += coefficients(idx++) * position(0);
        result += coefficients(idx++) * position(1);
        result += coefficients(idx++) * position(2);
    }
    if constexpr (F >= BILINEAR) {
        result += coefficients(idx++) * position(1) * position(2);
        result += coefficients(idx++) * position(0) * position(2);
        result += coefficients(idx++) * position(0) * position(1);
    }
    if constexpr (F >= QUADRATIC) {
        result += coefficients(idx++) * position(0) * position(0);
        result += coefficients(idx++) * position(1) * position(1);
        result += coefficients(idx++) * position(2) * position(2);
    }
    if constexpr (F >= BILINQUAD) {
        result += coefficients(idx++) * position(0) * position(1) * position(1);
        result += coefficients(idx++) * position(0) * position(2) * position(2);
        result += coefficients(idx++) * position(1) * position(0) * position(0);
        result += coefficients(idx++) * position(1) * position(2) * position(2);
        result += coefficients(idx++) * position(2) * position(0) * position(0);
        result += coefficients(idx++) * position(2) * position(1) * position(1);
        result += coefficients(idx++) * position(0) * position(1) * position(2);
    }
    if constexpr (F >= CUBIC) {
        result += coefficients(idx++) * position(0) * position(0) * position(0);
        result += coefficients(idx++) * position(1) * position(1) * position(1);
        result += coefficients(idx++) * position(2) * position(2) * position(2);
    }

    return result;
}

/**
 * @brief Computes predicted values for all nodes using the fitted coefficients.
 */
template<InterpolationFunction F>
DynamicMatrix compute_predicted_values(const DynamicMatrix& coefficients, const NodeData& xyz) {
    DynamicMatrix predicted(xyz.rows(), coefficients.cols());
    for (Index column = 0; column < coefficients.cols(); ++column) {
        for (Index row = 0; row < xyz.rows(); ++row) {
            predicted(row, column) = evaluate<F>(coefficients.col(column), xyz.row(row));
        }
    }
    return predicted;
}

/**
 * @brief Computes column-wise coefficients of determination (R²).
 */
DynamicVector compute_r2(const DynamicMatrix& predicted, const NodeData& values) {
    DynamicVector r2(values.cols());
    for (Index column = 0; column < values.cols(); ++column) {
        const Precision mean = values.col(column).mean();
        const Precision ss_tot = (values.col(column).array() - mean).square().sum();
        const Precision ss_res = (values.col(column).array() - predicted.col(column).array()).square().sum();
        r2(column) = ss_tot < 1e-10 ? Precision(1) : Precision(1) - ss_res / ss_tot;
    }
    return r2;
}

} // namespace

/**
 * @brief Template implementation of the polynomial interpolation routine.
 */
template<InterpolationFunction F>
DynamicMatrix interpolate(const NodeData& xyz,
                          const NodeData& values,
                          const Vec3& center,
                          DynamicVector* r2_values) {
    NodeData adjusted_xyz = xyz;
    for (Index row = 0; row < adjusted_xyz.rows(); ++row) {
        adjusted_xyz.row(row) -= center.transpose();
    }

    Vec3 adjusted_center = Vec3::Zero();

    SemiStaticMatrix<term_count<F>()> lhs{adjusted_xyz.rows(), term_count<F>()};
    for (Index row = 0; row < adjusted_xyz.rows(); ++row) {
        if constexpr (F >= CONSTANT) {
            lhs(row, 0) = 1.0;
        }
        if constexpr (F >= LINEAR) {
            lhs(row, 1) = adjusted_xyz(row, 0);
            lhs(row, 2) = adjusted_xyz(row, 1);
            lhs(row, 3) = adjusted_xyz(row, 2);
        }
        if constexpr (F >= BILINEAR) {
            lhs(row, 4) = adjusted_xyz(row, 1) * adjusted_xyz(row, 2);
            lhs(row, 5) = adjusted_xyz(row, 0) * adjusted_xyz(row, 2);
            lhs(row, 6) = adjusted_xyz(row, 0) * adjusted_xyz(row, 1);
        }
        if constexpr (F >= QUADRATIC) {
            lhs(row, 7) = adjusted_xyz(row, 0) * adjusted_xyz(row, 0);
            lhs(row, 8) = adjusted_xyz(row, 1) * adjusted_xyz(row, 1);
            lhs(row, 9) = adjusted_xyz(row, 2) * adjusted_xyz(row, 2);
        }
        if constexpr (F >= BILINQUAD) {
            lhs(row, 10) = adjusted_xyz(row, 0) * adjusted_xyz(row, 1) * adjusted_xyz(row, 1);
            lhs(row, 11) = adjusted_xyz(row, 0) * adjusted_xyz(row, 2) * adjusted_xyz(row, 2);
            lhs(row, 12) = adjusted_xyz(row, 1) * adjusted_xyz(row, 0) * adjusted_xyz(row, 0);
            lhs(row, 13) = adjusted_xyz(row, 1) * adjusted_xyz(row, 2) * adjusted_xyz(row, 2);
            lhs(row, 14) = adjusted_xyz(row, 2) * adjusted_xyz(row, 0) * adjusted_xyz(row, 0);
            lhs(row, 15) = adjusted_xyz(row, 2) * adjusted_xyz(row, 1) * adjusted_xyz(row, 1);
            lhs(row, 16) = adjusted_xyz(row, 0) * adjusted_xyz(row, 1) * adjusted_xyz(row, 2);
        }
        if constexpr (F >= CUBIC) {
            lhs(row, 17) = adjusted_xyz(row, 0) * adjusted_xyz(row, 0) * adjusted_xyz(row, 0);
            lhs(row, 18) = adjusted_xyz(row, 1) * adjusted_xyz(row, 1) * adjusted_xyz(row, 1);
            lhs(row, 19) = adjusted_xyz(row, 2) * adjusted_xyz(row, 2) * adjusted_xyz(row, 2);
        }
    }

    DynamicMatrix ata = lhs.transpose() * lhs;
    const Precision determinant = ata.determinant();

    if (determinant < 1e-10 || determinant > 1e20) {
        if (r2_values) {
            *r2_values = DynamicVector::Zero(values.cols());
        }
        return interpolate<demote_order<F>()>(xyz, values, center, r2_values);
    }

    DynamicMatrix atb = lhs.transpose() * values;
    DynamicMatrix coefficients{ata.rows(), atb.cols()};
    auto solver = ata.fullPivHouseholderQr();
    for (Index column = 0; column < values.cols(); ++column) {
        coefficients.col(column) = solver.solve(atb.col(column));
    }

    DynamicMatrix results(1, values.cols());
    for (Index column = 0; column < values.cols(); ++column) {
        results(0, column) = evaluate<F>(coefficients.col(column), adjusted_center);
    }

    for (Index column = 0; column < results.cols(); ++column) {
        for (Index row = 0; row < results.rows(); ++row) {
            if (std::isnan(results(row, column)) || std::isinf(results(row, column))) {
                return interpolate<demote_order<F>()>(xyz, values, center, r2_values);
            }
        }
    }

    if (r2_values) {
        const DynamicMatrix predicted = compute_predicted_values<F>(coefficients, adjusted_xyz);
        *r2_values = compute_r2(predicted, values);
    }

    return results;
}

#define INSTANTIATE(FUNC) \
    template DynamicMatrix interpolate<FUNC>(const NodeData&, const NodeData&, const Vec3&, DynamicVector*);

INSTANTIATE(CONSTANT)
INSTANTIATE(LINEAR)
INSTANTIATE(BILINEAR)
INSTANTIATE(QUADRATIC)
INSTANTIATE(BILINQUAD)
INSTANTIATE(CUBIC)

#undef INSTANTIATE

DynamicMatrix interpolate(const NodeData& xyz,
                          const NodeData& values,
                          const Vec3& center,
                          DynamicVector* r2_values,
                          float accuracy_factor,
                          InterpolationFunction max_accuracy) {
    const float weighted_rows = static_cast<float>(xyz.rows()) * accuracy_factor;

    if (weighted_rows < term_count<CONSTANT>() || max_accuracy == InterpolationFunction::CONSTANT) {
        return interpolate<CONSTANT>(xyz, values, center, r2_values);
    }
    if (weighted_rows < term_count<LINEAR>() || max_accuracy == InterpolationFunction::LINEAR) {
        return interpolate<LINEAR>(xyz, values, center, r2_values);
    }
    if (weighted_rows < term_count<BILINEAR>() || max_accuracy == InterpolationFunction::BILINEAR) {
        return interpolate<BILINEAR>(xyz, values, center, r2_values);
    }
    if (weighted_rows < term_count<QUADRATIC>() || max_accuracy == InterpolationFunction::QUADRATIC) {
        return interpolate<QUADRATIC>(xyz, values, center, r2_values);
    }
    if (weighted_rows < term_count<BILINQUAD>() || max_accuracy == InterpolationFunction::BILINQUAD) {
        return interpolate<BILINQUAD>(xyz, values, center, r2_values);
    }
    return interpolate<CUBIC>(xyz, values, center, r2_values);
}

Interpolator::Interpolator(InterpolationFunction method, float accuracy)
    : m_method(method), m_accuracy(accuracy) {}

void Interpolator::set_function(InterpolationFunction method) {
    m_method = method;
}

InterpolationFunction Interpolator::get_function() const {
    return m_method;
}

void Interpolator::set_accuracy(float accuracy) {
    m_accuracy = accuracy;
}

float Interpolator::get_accuracy() const {
    return m_accuracy;
}

DynamicMatrix Interpolator::operator()(const NodeData& xyz,
                                       const NodeData& values,
                                       const Vec3& center,
                                       DynamicVector* r2_values) {
    return interpolate(xyz, values, center, r2_values, m_accuracy, m_method);
}

} // namespace interpolate
} // namespace math
} // namespace fem

