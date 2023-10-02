//
// Created by Luecx on 30.09.2023.
//


#include "interpolate.h"
#include <iostream>

namespace fem {
namespace math {
namespace interpolate {

/**
 * \brief Determine the number of terms for a given interpolation function.
 *
 * \tparam F The interpolation function.
 * \return Number of terms for interpolation function F.
 */
template<InterpolationFunction F>
constexpr inline int n_terms() {
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
        return 0;    // Default case, should not be reached.
    }
}

template<InterpolationFunction F>
constexpr InterpolationFunction get_next_lower() {
    switch (F) {
        case CUBIC: return BILINQUAD;
        case BILINQUAD: return QUADRATIC;
        case QUADRATIC: return BILINEAR;
        case BILINEAR: return LINEAR;
        case LINEAR: return CONSTANT;
        case CONSTANT: return CONSTANT; // or throw an exception, etc.
        default: return CONSTANT;
    }
}

/**
 * \brief Evaluates an interpolation polynomial at a given point.
 *
 * \tparam F The interpolation function.
 * \param[in] coeff The coefficients of the interpolation polynomial.
 * \param[in] center The point where the polynomial is to be evaluated.
 * \return Resultant value after interpolation evaluation.
 *
 * \example
 * DynamicVector coefficients = {2, 3, 4};
 * StaticVector<3> point = {1, 2, 3};
 * Precision result = evaluate<LINEAR>(coefficients, point);
 */
template<InterpolationFunction F>
inline Precision evaluate(const DynamicVector& coeff, const StaticVector<3>& center) {
    Precision result = 0;

    int       idx    = 0;

    if constexpr (F >= CONSTANT) {
        result += coeff(idx++);
    }
    if constexpr (F >= LINEAR) {
        result += coeff(idx++) * center(0);
        result += coeff(idx++) * center(1);
        result += coeff(idx++) * center(2);
    }
    if constexpr (F >= BILINEAR) {
        result += coeff(idx++) * center(1) * center(2);
        result += coeff(idx++) * center(0) * center(2);
        result += coeff(idx++) * center(0) * center(1);
    }
    if constexpr (F >= QUADRATIC) {
        result += coeff(idx++) * center(0) * center(0);
        result += coeff(idx++) * center(1) * center(1);
        result += coeff(idx++) * center(2) * center(2);
    }
    if constexpr (F >= BILINQUAD) {
        result += coeff(idx++) * center(0) * center(1) * center(1);
        result += coeff(idx++) * center(0) * center(2) * center(2);

        result += coeff(idx++) * center(1) * center(0) * center(0);
        result += coeff(idx++) * center(1) * center(2) * center(2);

        result += coeff(idx++) * center(2) * center(0) * center(0);
        result += coeff(idx++) * center(2) * center(1) * center(1);

        result += coeff(idx++) * center(0) * center(1) * center(2);
    }
    if constexpr (F >= CUBIC) {
        result += coeff(idx++) * center(0) * center(0) * center(0);
        result += coeff(idx++) * center(1) * center(1) * center(1);
        result += coeff(idx++) * center(2) * center(2) * center(2);
    }

    return result;
}


/**
 * \brief Computes the predicted values for given coefficients and interpolation nodes.
 *
 * \tparam F The interpolation function.
 * \param[in] coe The coefficients matrix.
 * \param[in] xyz The interpolation nodes.
 * \return Matrix containing the predicted values.
 */
template<InterpolationFunction F>
DynamicMatrix compute_predicted_values(const DynamicMatrix& coe, const NodeData& xyz) {
    DynamicMatrix predicted_values(xyz.rows(), coe.cols());
    for (int i = 0; i < coe.cols(); i++) {
        for (int j = 0; j < xyz.rows(); j++) {
            predicted_values(j, i) = evaluate<F>(coe.col(i), xyz.row(j));
        }
    }
    return predicted_values;
}

/**
 * \brief Computes the coefficient of determination (R^2) for predicted vs actual values.
 *
 * \param[in] predicted_values The predicted values from interpolation.
 * \param[in] values The actual values.
 * \return Vector of R^2 values for each column of input matrices.
 */
DynamicVector compute_r2(const DynamicMatrix& predicted_values, const NodeData& values) {
    DynamicVector r2_values(values.cols());
    for (int i = 0; i < values.cols(); i++) {
        Precision mean   = values.col(i).mean();
        Precision ss_tot = (values.col(i).array() - mean).square().sum();
        Precision ss_res = (values.col(i).array() - predicted_values.col(i).array()).square().sum();

        // Check for very small variance
        if (ss_tot < 1e-10) {
            r2_values(i) = 1.0;    // If variance is negligible, assign R^2 to be 1.
        } else {
            r2_values(i) = 1 - ss_res / ss_tot;
        }
    }
    return r2_values;
}

template<InterpolationFunction F>
DynamicMatrix
    interpolate(const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, DynamicVector* r2_values) {
    
    NodeData adjusted_xyz = xyz;
    // Subtract center from each row of xyz
    for (int i = 0; i < xyz.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            adjusted_xyz(i, j) -= center(j);
        }
    }

    // Subtract center from center (which makes it a zero vector)
    StaticVector<3> adjusted_center = center;
    for (int j = 0; j < 3; ++j) {
        adjusted_center(j) -= center(j);
    }
    
    SemiStaticMatrix<n_terms<F>()> lhs {adjusted_xyz.rows(), n_terms<F>()};
    for (int r = 0; r < adjusted_xyz.rows(); r++) {
        if constexpr (F >= CONSTANT) {
            lhs(r, 0) = 1;
        }
        if constexpr (F >= LINEAR) {
            lhs(r, 1) = adjusted_xyz(r, 0);
            lhs(r, 2) = adjusted_xyz(r, 1);
            lhs(r, 3) = adjusted_xyz(r, 2);
        }
        if constexpr (F >= BILINEAR) {
            lhs(r, 4) = adjusted_xyz(r, 1) * adjusted_xyz(r, 2);
            lhs(r, 5) = adjusted_xyz(r, 0) * adjusted_xyz(r, 2);
            lhs(r, 6) = adjusted_xyz(r, 0) * adjusted_xyz(r, 1);
        }
        if constexpr (F >= QUADRATIC) {
            lhs(r, 7) = adjusted_xyz(r, 0) * adjusted_xyz(r, 0);
            lhs(r, 8) = adjusted_xyz(r, 1) * adjusted_xyz(r, 1);
            lhs(r, 9) = adjusted_xyz(r, 2) * adjusted_xyz(r, 2);
        }
        if constexpr (F >= BILINQUAD) {
            lhs(r, 10) = adjusted_xyz(r, 0) * adjusted_xyz(r, 1) * adjusted_xyz(r, 1);
            lhs(r, 11) = adjusted_xyz(r, 0) * adjusted_xyz(r, 2) * adjusted_xyz(r, 2);

            lhs(r, 12) = adjusted_xyz(r, 1) * adjusted_xyz(r, 0) * adjusted_xyz(r, 0);
            lhs(r, 13) = adjusted_xyz(r, 1) * adjusted_xyz(r, 2) * adjusted_xyz(r, 2);

            lhs(r, 14) = adjusted_xyz(r, 2) * adjusted_xyz(r, 0) * adjusted_xyz(r, 0);
            lhs(r, 15) = adjusted_xyz(r, 2) * adjusted_xyz(r, 1) * adjusted_xyz(r, 1);

            lhs(r, 16) = adjusted_xyz(r, 0) * adjusted_xyz(r, 1) * adjusted_xyz(r, 2);
        }
        if constexpr (F >= CUBIC) {
            lhs(r, 17) = adjusted_xyz(r, 0) * adjusted_xyz(r, 0) * adjusted_xyz(r, 0);
            lhs(r, 18) = adjusted_xyz(r, 1) * adjusted_xyz(r, 1) * adjusted_xyz(r, 1);
            lhs(r, 19) = adjusted_xyz(r, 2) * adjusted_xyz(r, 2) * adjusted_xyz(r, 2);
        }
    }
    DynamicMatrix results(1, values.cols());

    auto ATA = lhs.transpose() * lhs;
    auto determinant = ATA.determinant();
    auto ATB = lhs.transpose() * values;
    auto coe = ATA.fullPivHouseholderQr().solve(ATB);

    for (int i = 0; i < values.cols(); i++) {
        results(0, i) = evaluate<F>(coe.col(i), adjusted_center);
    }

    if (r2_values) {
        DynamicMatrix predicted_vals = compute_predicted_values<F>(coe, adjusted_xyz);
        *r2_values                   = compute_r2(predicted_vals, values);
    }

    return results;
}

template DynamicMatrix interpolate<CONSTANT> (const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, DynamicVector* r2_values);
template DynamicMatrix interpolate<LINEAR>   (const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, DynamicVector* r2_values);
template DynamicMatrix interpolate<BILINEAR> (const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, DynamicVector* r2_values);
template DynamicMatrix interpolate<QUADRATIC>(const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, DynamicVector* r2_values);
template DynamicMatrix interpolate<BILINQUAD>(const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, DynamicVector* r2_values);
template DynamicMatrix interpolate<CUBIC>    (const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, DynamicVector* r2_values);

DynamicMatrix interpolate(const NodeData&        xyz,
                          const NodeData&        values,
                          const StaticVector<3>& center,
                          DynamicVector*         r2_values,
                          float                  accuracy_factor,
                          InterpolationFunction  max_accuracy) {
    float adjusted_rows = xyz.rows() * accuracy_factor;

    if (adjusted_rows < n_terms<CONSTANT>() || max_accuracy == InterpolationFunction::CONSTANT) {
        return interpolate<CONSTANT>(xyz, values, center, r2_values);
    } else if (adjusted_rows < n_terms<LINEAR>() || max_accuracy == InterpolationFunction::LINEAR) {
        return interpolate<LINEAR>(xyz, values, center, r2_values);
    } else if (adjusted_rows < n_terms<BILINEAR>() || max_accuracy == InterpolationFunction::BILINEAR) {
        return interpolate<BILINEAR>(xyz, values, center, r2_values);
    } else if (adjusted_rows < n_terms<QUADRATIC>() || max_accuracy == InterpolationFunction::QUADRATIC) {
        return interpolate<QUADRATIC>(xyz, values, center, r2_values);
    } else if (adjusted_rows < n_terms<BILINQUAD>() || max_accuracy == InterpolationFunction::BILINQUAD) {
        return interpolate<BILINQUAD>(xyz, values, center, r2_values);
    } else {
        return interpolate<CUBIC>(xyz, values, center, r2_values);
    }
}

Interpolator::Interpolator(InterpolationFunction method_, float accuracy_)
    : method(method_), accuracy(accuracy_) {}
void Interpolator::set_function(InterpolationFunction method_) {
    method = method_;
}
InterpolationFunction Interpolator::get_function() const {
    return method;
}
void Interpolator::set_accuracy(float accuracy_) {
    accuracy = accuracy_;
}
float Interpolator::get_accuracy() const {
    return accuracy;
}
DynamicMatrix Interpolator::operator()(const NodeData&        xyz,
                                       const NodeData&        values,
                                       const StaticVector<3>& center,
                                       DynamicVector*         r2_values) {
    return interpolate(xyz, values, center, r2_values, accuracy, method);
}
}    // namespace interpolate
}    // namespace math
}    // namespace fem