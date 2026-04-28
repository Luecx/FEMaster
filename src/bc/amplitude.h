/**
 * @file amplitude.h
 * @brief Declares time-dependent amplitude series for boundary conditions.
 *
 * Amplitudes provide reusable scalar time histories that can be assigned to
 * loads or other boundary conditions. They evaluate to a scalar multiplier at a
 * requested time based on the configured interpolation method.
 *
 * @see src/bc/amplitude.cpp
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "../core/types_num.h"
#include "../data/namable.h"

#include <memory>
#include <string>
#include <vector>

namespace fem {
namespace bc {
/**
 * @enum Interpolation
 * @brief Enumerates interpolation schemes for amplitude evaluation.
 */
enum class Interpolation {
    Step,    ///< Holds the previous sample until the next sample is reached.
    Nearest, ///< Selects the sample closest to the query time.
    Linear   ///< Linear interpolation between adjacent samples.
};

/**
 * @struct Amplitude
 * @brief Stores a time-value series with configurable interpolation.
 */
struct Amplitude : Namable {
    using Ptr = std::shared_ptr<Amplitude>; ///< Shared pointer alias for storage.

    /**
     * @struct Sample
     * @brief Stores one amplitude support point.
     */
    struct Sample {
        Precision time_  = 0; ///< Sample time coordinate.
        Precision value_ = 0; ///< Sample value at the given time.
    };

    /**
     * @brief Creates an amplitude with a name and interpolation mode.
     *
     * @param name Amplitude identifier.
     * @param interpolation Interpolation used between stored samples.
     */
    explicit Amplitude(const std::string& name = "",
                       Interpolation interpolation = Interpolation::Linear);

    /**
     * @brief Changes how values between samples are evaluated.
     *
     * @param interpolation New interpolation mode.
     */
    void set_interpolation(Interpolation interpolation);

    /**
     * @brief Returns the active interpolation mode.
     *
     * @return Interpolation Current interpolation mode.
     */
    [[nodiscard]] Interpolation interpolation() const;

    /**
     * @brief Removes all stored samples.
     *
     * Used when an existing named amplitude is redefined.
     */
    void clear_samples();

    /**
     * @brief Adds or replaces a time-value sample.
     *
     * Samples are kept sorted by time. A new sample with an already existing
     * time replaces the old value.
     *
     * @param time Sample time coordinate.
     * @param value Sample value.
     */
    void add_sample(Precision time, Precision value);

    /**
     * @brief Evaluates the amplitude at the requested time.
     *
     * Empty amplitudes evaluate to `1.0`, so a missing scale history behaves
     * like a neutral load multiplier.
     *
     * @param time Query time.
     * @return Precision Interpolated scalar multiplier.
     */
    [[nodiscard]] Precision evaluate(Precision time) const;

private:
    Interpolation       interpolation_;
    std::vector<Sample> samples_;
};
} // namespace bc
} // namespace fem
