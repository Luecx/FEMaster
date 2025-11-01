/**
 * @file amplitude.cpp
 * @brief Implements the time-dependent amplitude series utilities.
 */

#include "amplitude.h"

#include <cmath>
#include <limits>

namespace fem {
namespace bc {

namespace {
constexpr Precision kTimeEqualityEps = 1e-12;
}

Amplitude::Amplitude(const std::string& name, Interpolation interpolation)
    : Namable(name), interpolation_(interpolation) {}

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
    auto it = std::lower_bound(samples_.begin(), samples_.end(), time,
        [](const Sample& sample, Precision t) { return sample.time < t; });

    if (it != samples_.end() && std::abs(it->time - time) <= kTimeEqualityEps) {
        it->value = value;
        return;
    }

    samples_.insert(it, Sample{time, value});
}

Precision Amplitude::evaluate(Precision time) const {
    if (samples_.empty()) {
        return 1.0;
    }

    if (samples_.size() == 1) {
        return samples_.front().value;
    }

    auto it = std::lower_bound(samples_.begin(), samples_.end(), time,
        [](const Sample& sample, Precision t) { return sample.time < t; });

    switch (interpolation_) {
        case Interpolation::Step: {
            if (it == samples_.end()) {
                return samples_.back().value;
            }
            if (std::abs(it->time - time) <= kTimeEqualityEps) {
                return it->value;
            }
            if (it == samples_.begin()) {
                return it->value;
            }
            return std::prev(it)->value;
        }
        case Interpolation::Nearest: {
            if (it == samples_.end()) {
                return samples_.back().value;
            }
            if (it == samples_.begin()) {
                return it->value;
            }
            if (std::abs(it->time - time) <= kTimeEqualityEps) {
                return it->value;
            }
            const auto prev = std::prev(it);
            const Precision dist_prev = std::abs(time - prev->time);
            const Precision dist_next = std::abs(it->time - time);
            return dist_prev <= dist_next ? prev->value : it->value;
        }
        case Interpolation::Linear:
        default: {
            if (it == samples_.end()) {
                return samples_.back().value;
            }
            if (std::abs(it->time - time) <= kTimeEqualityEps || it == samples_.begin()) {
                return it->value;
            }
            const auto prev = std::prev(it);
            const Precision span = it->time - prev->time;
            if (std::abs(span) <= std::numeric_limits<Precision>::epsilon()) {
                return it->value;
            }
            const Precision alpha = (time - prev->time) / span;
            return prev->value + alpha * (it->value - prev->value);
        }
    }
}

} // namespace bc
} // namespace fem

