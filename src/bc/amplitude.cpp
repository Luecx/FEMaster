/**
 * @file amplitude.cpp
 * @brief Implements sorted amplitude storage and scalar time interpolation.
 *
 * This file contains the complete evaluation path used by `Amplitude`: tolerant
 * comparison of time coordinates, binary search in the sorted sample vector,
 * endpoint clamping and the three supported interpolation rules. It also
 * implements insertion with replacement of coincident support points so every
 * amplitude remains ordered and free of duplicate times.
 *
 * @see amplitude.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "amplitude.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {
namespace bc {
namespace {

// Absolute tolerance used only when deciding whether two stored or queried time
// coordinates represent the same nominal support point. The value avoids
// duplicate samples caused by insignificant parser or floating-point noise.
constexpr Precision kTimeEqualityEps = 1e-12;

/**
 * Compares two amplitude support-point times using the module-wide absolute
 * tolerance.
 *
 * Treating nearly equal coordinates as identical prevents duplicate support
 * points caused by parser conversion or floating-point round-off.
 *
 * @param lhs First time coordinate.
 * @param rhs Second time coordinate.
 * @return True when the absolute difference is within the configured tolerance.
 */
bool same_time(Precision lhs, Precision rhs) {
    return std::abs(lhs - rhs) <= kTimeEqualityEps;
}

/**
 * Locates the first amplitude sample at or after a requested analysis time.
 *
 * Samples are maintained in ascending time order, so binary search finds the
 * upper interpolation bracket in logarithmic time. The returned iterator may
 * equal samples.end() when the query is beyond the final support point.
 *
 * @param samples Sorted amplitude support points.
 * @param time Requested analysis time.
 * @return Iterator to the first sample whose time is not less than time.
 */
std::vector<Amplitude::Sample>::const_iterator lower_sample(const std::vector<Amplitude::Sample>& samples,
                                                            Precision time) {
    return std::lower_bound(samples.begin(), samples.end(), time,
        [](const Amplitude::Sample& sample, Precision sample_time) {
            // The comparator orders a stored sample against the scalar search
            // key using only its time coordinate.
            return sample.time_ < sample_time;
        });
}

/**
 * Evaluates a sampled amplitude using left-continuous step interpolation.
 *
 * Queries before and after the sampled interval are clamped to the endpoint
 * values. Interior queries retain the preceding value until the next support
 * time is reached.
 *
 * @param samples Sorted amplitude support points.
 * @param time Requested analysis time.
 * @return The step-interpolated amplitude value.
 */
Precision evaluate_step(const std::vector<Amplitude::Sample>& samples, Precision time) {
    // `upper` is either the first sample at or after the query or `end()` when
    // the query lies beyond the complete sampled interval.
    const auto upper = lower_sample(samples, time);

    // Hold the final value for all times after the last support point.
    if (upper == samples.end()) {
        return samples.back().value_;
    }

    // Clamp times before the first sample to its value. At an exact support
    // point, activate the value stored at that support point immediately.
    if (same_time(upper->time_, time) || upper == samples.begin()) {
        return upper->value_;
    }

    // For a strict interior query, the previous sample supplies the held value.
    return std::prev(upper)->value_;
}

/**
 * Evaluates a sampled amplitude using nearest-neighbor interpolation.
 *
 * Endpoint values are used outside the sampled interval. Within the interval,
 * the two bracketing samples are compared by temporal distance.
 *
 * @param samples Sorted amplitude support points.
 * @param time Requested analysis time.
 * @return The value of the nearest support point.
 */
Precision evaluate_nearest(const std::vector<Amplitude::Sample>& samples, Precision time) {
    // Find the first sample not earlier than the requested time.
    const auto upper = lower_sample(samples, time);

    // Queries beyond the final support point use the final stored value.
    if (upper == samples.end()) {
        return samples.back().value_;
    }

    // Queries before the first point and queries coincident with a point require
    // no distance comparison.
    if (upper == samples.begin() || same_time(upper->time_, time)) {
        return upper->value_;
    }

    // Compare the absolute distances to the lower and upper neighbors. The
    // `<=` relation intentionally resolves an exact midpoint toward the earlier
    // sample, producing deterministic behavior.
    const auto      lower     = std::prev(upper);
    const Precision dist_low  = std::abs(time - lower->time_);
    const Precision dist_high = std::abs(upper->time_ - time);
    return dist_low <= dist_high ? lower->value_ : upper->value_;
}

/**
 * Evaluates a sampled amplitude by affine interpolation between support points.
 *
 * Values outside the sampled interval are clamped to the nearest endpoint.
 * Interior queries use the normalized time coordinate between the two
 * bracketing samples.
 *
 * @param samples Sorted amplitude support points.
 * @param time Requested analysis time.
 * @return The linearly interpolated amplitude value.
 */
Precision evaluate_linear(const std::vector<Amplitude::Sample>& samples, Precision time) {
    // Find the upper bracket of the interpolation interval.
    const auto upper = lower_sample(samples, time);

    // Clamp queries after the complete history to the final value.
    if (upper == samples.end()) {
        return samples.back().value_;
    }

    // Clamp queries before the first sample and return exact support-point
    // values directly to avoid unnecessary arithmetic.
    if (upper == samples.begin() || same_time(upper->time_, time)) {
        return upper->value_;
    }

    // The previous iterator and `upper` delimit the interpolation interval.
    const auto      lower = std::prev(upper);
    const Precision span  = upper->time_ - lower->time_;

    // Duplicate times should already be prevented by `add_sample()`. The guard
    // nevertheless avoids division by a numerically zero interval if malformed
    // data reaches this helper.
    if (std::abs(span) <= std::numeric_limits<Precision>::epsilon()) {
        return upper->value_;
    }

    // Normalize the query into the interval and interpolate the sample values
    // with the standard affine combination.
    const Precision alpha = (time - lower->time_) / span;
    return lower->value_ + alpha * (upper->value_ - lower->value_);
}
} // namespace

/**
 * Constructs an amplitude definition with a name and interpolation rule.
 *
 * The sample history starts empty, making the amplitude neutral until support
 * points are inserted. The selected rule is retained for later evaluations.
 *
 * @param name Name used in model data and diagnostics.
 * @param interpolation Rule used between stored support points.
 */
Amplitude::Amplitude(const std::string& name, Interpolation interpolation)
    : Namable(name),
      interpolation_(interpolation) {}

/**
 * Selects the interpolation rule used for future amplitude evaluations.
 *
 * Changing the rule preserves the existing, sorted support-point history and
 * only changes how values between those points are computed.
 *
 * @param interpolation New interpolation rule.
 */
void Amplitude::set_interpolation(Interpolation interpolation) {
    interpolation_ = interpolation;
}

/**
 * Returns the interpolation rule currently associated with this amplitude.
 *
 * This is the rule used by evaluate for histories containing multiple support
 * points.
 *
 * @return The configured interpolation rule.
 */
Interpolation Amplitude::interpolation() const {
    return interpolation_;
}

/**
 * Removes all support points from the amplitude history.
 *
 * The name and interpolation rule remain unchanged. Subsequent evaluations use
 * the neutral multiplier associated with an empty optional amplitude.
 */
void Amplitude::clear_samples() {
    samples_.clear();
}

/**
 * Inserts or replaces one amplitude support point while preserving time order.
 *
 * A sample within the configured equality tolerance is replaced; otherwise a
 * new sample is inserted at the binary-search position.
 *
 * @param time Support-point time coordinate.
 * @param value Amplitude multiplier at the support point.
 */
void Amplitude::add_sample(Precision time, Precision value) {
    // Find the first sample whose time is not smaller than the new coordinate.
    // This iterator is both the replacement candidate and the correct insertion
    // position when no coincident sample exists.
    auto insert_position = std::lower_bound(samples_.begin(), samples_.end(), time,
        [](const Sample& sample, Precision sample_time) {
            return sample.time_ < sample_time;
        });

    // A repeated definition replaces the stored value in place. Keeping only
    // one sample per nominal time prevents zero-length interpolation intervals.
    if (insert_position != samples_.end() && same_time(insert_position->time_, time)) {
        insert_position->value_ = value;
        return;
    }

    // Insert before the first later sample so the vector remains sorted after
    // every update.
    samples_.insert(insert_position, Sample{time, value});
}

/**
 * Evaluates the amplitude history at a requested analysis time.
 *
 * Empty and single-point histories are handled as constant histories. Larger
 * histories dispatch to the configured step, nearest-neighbor or linear rule.
 *
 * @param time Requested analysis time.
 * @return The amplitude multiplier applied to the associated load.
 */
Precision Amplitude::evaluate(Precision time) const {
    // An empty optional amplitude must leave the nominal load unchanged.
    if (samples_.empty()) {
        return 1.0;
    }

    // A single support point defines a constant history for all times.
    if (samples_.size() == 1) {
        return samples_.front().value_;
    }

    // Dispatch to the interpolation-specific helper. Every helper receives a
    // sorted vector with at least two entries and performs endpoint clamping.
    switch (interpolation_) {
        case Interpolation::Step:
            return evaluate_step(samples_, time);
        case Interpolation::Nearest:
            return evaluate_nearest(samples_, time);
        case Interpolation::Linear:
            return evaluate_linear(samples_, time);
    }

    // Preserve a deterministic result if an invalid enum value is introduced
    // through external data or a future extension without a matching case.
    return evaluate_linear(samples_, time);
}
} // namespace bc
} // namespace fem
