/**
 * @file amplitude.cpp
 * @brief Implements the time-dependent amplitude series utilities.
 */

#include "amplitude.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {
namespace bc {
namespace {
constexpr Precision kTimeEqualityEps = 1e-12;

/**
 * @brief Returns true when two time coordinates should be treated as equal.
 *
 * Amplitudes are usually parsed from text input. A small absolute tolerance
 * avoids creating duplicate support points when the same nominal time appears
 * with tiny floating-point noise.
 *
 * @param lhs First time value.
 * @param rhs Second time value.
 * @return bool True when the values are within the time equality tolerance.
 */
bool same_time(Precision lhs, Precision rhs) {
    return std::abs(lhs - rhs) <= kTimeEqualityEps;
}

/**
 * @brief Finds the first sample whose time is not less than the query time.
 *
 * The sample list is always sorted by `add_sample`, so all interpolation modes
 * can use this one lookup and then inspect the neighboring sample if needed.
 *
 * @param samples Sorted sample list.
 * @param time Query time.
 * @return std::vector<Amplitude::Sample>::const_iterator Lower-bound iterator.
 */
std::vector<Amplitude::Sample>::const_iterator lower_sample(const std::vector<Amplitude::Sample>& samples,
                                                            Precision time) {
    return std::lower_bound(samples.begin(), samples.end(), time,
                            [](const Amplitude::Sample& sample, Precision sample_time) {
                                return sample.time_ < sample_time;
                            });
}

/**
 * @brief Evaluates left-continuous step interpolation.
 *
 * Before the first sample the first sample value is used. After the final sample
 * the final value is held. Between samples the previous sample is held.
 *
 * @param samples Sorted sample list with at least two entries.
 * @param time Query time.
 * @return Precision Interpolated value.
 */
Precision evaluate_step(const std::vector<Amplitude::Sample>& samples, Precision time) {
    const auto upper = lower_sample(samples, time);

    if (upper == samples.end()) {
        return samples.back().value_;
    }
    if (same_time(upper->time_, time) || upper == samples.begin()) {
        return upper->value_;
    }
    return std::prev(upper)->value_;
}

/**
 * @brief Evaluates nearest-neighbor interpolation.
 *
 * Boundary queries clamp to the first or last sample. Interior queries compare
 * the distance to the previous and next sample; ties prefer the previous value
 * for deterministic behavior.
 *
 * @param samples Sorted sample list with at least two entries.
 * @param time Query time.
 * @return Precision Interpolated value.
 */
Precision evaluate_nearest(const std::vector<Amplitude::Sample>& samples, Precision time) {
    const auto upper = lower_sample(samples, time);

    if (upper == samples.end()) {
        return samples.back().value_;
    }
    if (upper == samples.begin() || same_time(upper->time_, time)) {
        return upper->value_;
    }

    const auto      lower     = std::prev(upper);
    const Precision dist_low  = std::abs(time - lower->time_);
    const Precision dist_high = std::abs(upper->time_ - time);
    return dist_low <= dist_high ? lower->value_ : upper->value_;
}

/**
 * @brief Evaluates linear interpolation.
 *
 * Outside the sampled range the nearest endpoint value is used. Inside the
 * range the function interpolates between the previous and next support point.
 *
 * @param samples Sorted sample list with at least two entries.
 * @param time Query time.
 * @return Precision Interpolated value.
 */
Precision evaluate_linear(const std::vector<Amplitude::Sample>& samples, Precision time) {
    const auto upper = lower_sample(samples, time);

    if (upper == samples.end()) {
        return samples.back().value_;
    }
    if (upper == samples.begin() || same_time(upper->time_, time)) {
        return upper->value_;
    }

    const auto      lower = std::prev(upper);
    const Precision span  = upper->time_ - lower->time_;
    if (std::abs(span) <= std::numeric_limits<Precision>::epsilon()) {
        return upper->value_;
    }

    const Precision alpha = (time - lower->time_) / span;
    return lower->value_ + alpha * (upper->value_ - lower->value_);
}
} // namespace

Amplitude::Amplitude(const std::string& name, Interpolation interpolation)
    : Namable(name),
      interpolation_(interpolation) {}

void Amplitude::set_interpolation(Interpolation interpolation) {
    interpolation_ = interpolation;
}

Interpolation Amplitude::interpolation() const {
    return interpolation_;
}

void Amplitude::clear_samples() {
    samples_.clear();
}

void Amplitude::add_sample(Precision time, Precision value) {
    auto insert_position = std::lower_bound(samples_.begin(), samples_.end(), time,
                                            [](const Sample& sample, Precision sample_time) {
                                                return sample.time_ < sample_time;
                                            });

    // Repeated definitions at the same time replace the old value. This lets a
    // parser redefine a point without leaving duplicate support points behind.
    if (insert_position != samples_.end() && same_time(insert_position->time_, time)) {
        insert_position->value_ = value;
        return;
    }

    samples_.insert(insert_position, Sample{time, value});
}

Precision Amplitude::evaluate(Precision time) const {
    // A missing amplitude table must not change the load. Returning one keeps
    // optional amplitude use transparent for all load types.
    if (samples_.empty()) {
        return 1.0;
    }

    if (samples_.size() == 1) {
        return samples_.front().value_;
    }

    switch (interpolation_) {
        case Interpolation::Step:
            return evaluate_step(samples_, time);
        case Interpolation::Nearest:
            return evaluate_nearest(samples_, time);
        case Interpolation::Linear:
            return evaluate_linear(samples_, time);
    }

    return evaluate_linear(samples_, time);
}
} // namespace bc
} // namespace fem
