/**
 * @file amplitude.h
 * @brief Declares reusable scalar amplitude histories for boundary conditions.
 *
 * An `Amplitude` stores a named, time-ordered sequence of scalar support points
 * and evaluates the corresponding multiplier at an arbitrary analysis time.
 * Loads can reference the same amplitude object to share a prescribed temporal
 * evolution without duplicating interpolation or storage logic. The concrete
 * interpolation algorithms, sorted insertion and endpoint handling are
 * implemented in `amplitude.cpp`.
 *
 * @see Amplitude
 * @see Interpolation
 * @see amplitude.cpp
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
 * @brief Selects how an amplitude is evaluated between stored support points.
 *
 * All interpolation modes clamp queries outside the sampled time range to the
 * nearest endpoint value. They differ only in how an interior query is mapped
 * to the neighboring samples.
 */
enum class Interpolation {
    // Hold the value of the preceding sample until the next sample time is
    // reached. At an exact sample time the new sample value becomes active.
    Step,

    // Select the temporally closest sample. Equal distances are resolved in
    // favor of the preceding sample to keep the result deterministic.
    Nearest,

    // Interpolate affinely between the samples directly below and above the
    // query time.
    Linear
};

/**
 * @brief Stores and evaluates a named scalar time history.
 *
 * Samples are maintained in ascending time order. Adding a sample at an
 * already represented time replaces the existing value instead of introducing
 * a duplicate support point. An empty amplitude evaluates to the neutral
 * multiplier `1.0`, allowing optional amplitudes to be used without special
 * handling in every load implementation.
 */
struct Amplitude : Namable {
    // Shared ownership type used by loads and input-data collectors. Multiple
    // boundary conditions may intentionally reference the same time history.
    using Ptr = std::shared_ptr<Amplitude>;

    /**
     * @brief Represents one scalar support point of the amplitude history.
     *
     * The pair is kept intentionally lightweight because amplitude evaluation
     * searches a contiguous vector of samples.
     */
    struct Sample {
        // Analysis-time coordinate at which the stored value is prescribed.
        Precision time_  = 0;

        // Scalar amplitude value associated with `time_`.
        Precision value_ = 0;
    };

    // Construct a named amplitude with the requested interpolation rule. The
    // sample sequence initially remains empty and therefore evaluates to 1.0.
    explicit Amplitude(const std::string& name = "",
                       Interpolation interpolation = Interpolation::Linear);

    // Replace the interpolation rule used by subsequent evaluations without
    // modifying or reordering the stored support points.
    void set_interpolation(Interpolation interpolation);

    // Return the interpolation rule currently selected for this amplitude.
    [[nodiscard]] Interpolation interpolation() const;

    // Remove all support points while preserving the amplitude name and the
    // selected interpolation rule. This is used when redefining named input.
    void clear_samples();

    // Insert a support point into the sorted sequence. A time that matches an
    // existing sample within the implementation tolerance updates that sample's
    // value instead of creating a second entry.
    void add_sample(Precision time, Precision value);

    // Evaluate the scalar multiplier at `time`. Queries outside the sampled
    // interval are clamped to the nearest endpoint; an empty history returns
    // the neutral value 1.0.
    [[nodiscard]] Precision evaluate(Precision time) const;

private:
    // Active interpolation rule used for interior time queries.
    Interpolation interpolation_;

    // Sorted support-point storage maintained by `add_sample()`.
    std::vector<Sample> samples_;
};
} // namespace bc
} // namespace fem
