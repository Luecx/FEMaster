#pragma once

#include "../core/types.h"

#include <functional>
#include <ostream>

namespace fem {
namespace math {
namespace interpolate {

enum InterpolationFunction{
    CONSTANT,
    LINEAR,
    BILINEAR,
    QUADRATIC,
    BILINQUAD,
    CUBIC,
};

template<InterpolationFunction F>
constexpr int n_terms(){
    if constexpr (F == CONSTANT) {
        return 1;
    } else if constexpr (F == LINEAR) {
        return 4;
    } else if constexpr (F == BILINEAR) {  // Trilinear in 3D context
        return 7;
    } else if constexpr (F == QUADRATIC) {
        return 10;
    } else if constexpr (F == BILINQUAD) {
        return 10;
    } else if constexpr (F == CUBIC) {
        return 20;
    } else {
        return 0;  // Default case, should not be reached.
    }
}

template<InterpolationFunction F>
Precision evaluate(const DynamicVector& coeff, const StaticVector<3>& center) {
    Precision result = 0;

    int idx = 0;

    if constexpr (F >= CONSTANT) {
        result += coeff(idx++);
    }
    if constexpr(F >= LINEAR){
        result += coeff(idx++) * center(0);
        result += coeff(idx++) * center(1);
        result += coeff(idx++) * center(2);
    }
    if constexpr(F >= BILINEAR){
        result += coeff(idx++) * center(1) * center(2);
        result += coeff(idx++) * center(0) * center(2);
        result += coeff(idx++) * center(0) * center(1);
    }
    if constexpr(F >= QUADRATIC){
        result += coeff(idx++) * center(0) * center(0);
        result += coeff(idx++) * center(1) * center(1);
        result += coeff(idx++) * center(2) * center(2);
    }
    if constexpr(F >= BILINQUAD) {
        result += coeff(idx++) * center(0) * center(1) * center(1);
        result += coeff(idx++) * center(0) * center(2) * center(2);

        result += coeff(idx++) * center(1) * center(0) * center(0);
        result += coeff(idx++) * center(1) * center(2) * center(2);

        result += coeff(idx++) * center(2) * center(0) * center(0);
        result += coeff(idx++) * center(2) * center(1) * center(1);

        result += coeff(idx++) * center(0) * center(1) * center(2);
    }
    if constexpr(F >= CUBIC){
        result += coeff(idx++) * center(0) * center(0) * center(0);
        result += coeff(idx++) * center(1) * center(1) * center(1);
        result += coeff(idx++) * center(2) * center(2) * center(2);
    }

    return result;
}


template<InterpolationFunction F>
DynamicVector interpolate(const NodeData& xyz, const NodeData& values, const StaticVector<3>& center, Precision* buffer){
    MapMatrix lhs{buffer, xyz.rows(), n_terms<F>()};
    for (int r = 0; r < xyz.rows(); r++){
        if constexpr (F >= CONSTANT) {
            lhs(r,0) = 1;
        }
        if constexpr(F >= LINEAR){
            lhs(r,1) = xyz(r,0);
            lhs(r,2) = xyz(r,1);
            lhs(r,3) = xyz(r,2);
        }
        if constexpr(F >= BILINEAR){
            lhs(r,4) = xyz(r,1) * xyz(r,2);
            lhs(r,5) = xyz(r,0) * xyz(r,2);
            lhs(r,6) = xyz(r,0) * xyz(r,1);
        }
        if constexpr(F >= QUADRATIC){
            lhs(r,7) = xyz(r,0) * xyz(r,0);
            lhs(r,8) = xyz(r,1) * xyz(r,1);
            lhs(r,9) = xyz(r,2) * xyz(r,2);
        }
        if constexpr(F >= BILINQUAD) {
            lhs(r,10) = xyz(r,0) * xyz(r,1) * xyz(r,1);
            lhs(r,11) = xyz(r,0) * xyz(r,2) * xyz(r,2);

            lhs(r,12) = xyz(r,1) * xyz(r,0) * xyz(r,0);
            lhs(r,13) = xyz(r,1) * xyz(r,2) * xyz(r,2);

            lhs(r,14) = xyz(r,2) * xyz(r,0) * xyz(r,0);
            lhs(r,15) = xyz(r,2) * xyz(r,1) * xyz(r,1);

            lhs(r,16) = xyz(r,0) * xyz(r,1) * xyz(r,2);
        }
        if constexpr(F >= CUBIC){
            lhs(r,17) = xyz(r,0) * xyz(r,0) * xyz(r,0);
            lhs(r,18) = xyz(r,1) * xyz(r,1) * xyz(r,1);
            lhs(r,19) = xyz(r,2) * xyz(r,2) * xyz(r,2);
        }
    }
    DynamicVector results(values.cols());
    for (int col = 0; col < values.cols(); ++col) {
        auto sol = lhs.colPivHouseholderQr().solve(values.col(col));
        results(col) = evaluate<F>(sol, center);
    }
    return results;
}


}    // namespace interpolate
}    // namespace math
}    // namespace fem