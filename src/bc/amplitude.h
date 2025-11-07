wsl/**
 * @file amplitude.h
 * @brief Declares time-dependent amplitude series for boundary conditions.
 *
 * Amplitudes provide reusable scalar time histories that can be assigned to
 * loads or other boundary conditions. They evaluate to a scalar multiplier at a
 * requested time based on the configured interpolation method.
 */

#pragma once

#include "../core/types_num.h"
#include "../data/namable.h"

#include <algorithm>
#include <memory>
#include <vector>

namespace fem {
namespace bc {

/**
 * @enum Interpolation
 * @brief Enumerates interpolation schemes for amplitude evaluation.
 */
enum class Interpolation {
    Step,    ///< Left-continuous piecewise constant interpolation.
    Nearest, ///< Selects the sample closest to the query time.
    Linear   ///< Linear interpolation between adjacent samples.
};

/**
 * @struct Amplitude
 * @brief Stores a time-value series with configurable interpolation.
 */
struct Amplitude : Namable {
    using Ptr = std::shared_ptr<Amplitude>; ///< Shared pointer alias for storage.

    struct Sample {
        Precision time{0}; ///< Sample time coordinate.
        Precision value{0}; ///< Sample value at the given time.
    };

    explicit Amplitude(const std::string& name = "",
                       Interpolation interpolation = Interpolation::Linear);

    void set_interpolation(Interpolation interpolation);
    [[nodiscard]] Interpolation interpolation() const;

    void clear_samples();
    void add_sample(Precision time, Precision value);

    [[nodiscard]] Precision evaluate(Precision time) const;

private:
    Interpolation interpolation_;
    std::vector<Sample> samples_;
};

} // namespace bc
} // namespace fem
