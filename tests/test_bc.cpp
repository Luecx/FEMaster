#include "../src/bc/amplitude.h"

#include <gtest/gtest.h>

using namespace fem;

// 23) Amplitude insert + interpolation modes
TEST(BC_Amplitude, SamplesAndInterpolation) {
    bc::Amplitude amp("A", bc::Interpolation::Linear);
    // empty defaults to 1
    EXPECT_NEAR(amp.evaluate(0.0), 1.0, 1e-12);

    amp.add_sample(0.0, 0.0);
    amp.add_sample(1.0, 1.0);
    // overwrite existing time
    amp.add_sample(1.0, 2.0);

    // Linear
    amp.set_interpolation(bc::Interpolation::Linear);
    EXPECT_NEAR(amp.evaluate(0.0), 0.0, 1e-12);
    EXPECT_NEAR(amp.evaluate(1.0), 2.0, 1e-12);
    EXPECT_NEAR(amp.evaluate(0.5), 1.0, 1e-12);

    // Step (left-continuous): hold previous value until next sample
    amp.set_interpolation(bc::Interpolation::Step);
    EXPECT_NEAR(amp.evaluate(0.0), 0.0, 1e-12);
    EXPECT_NEAR(amp.evaluate(0.499999), 0.0, 1e-12);
    EXPECT_NEAR(amp.evaluate(0.5), 0.0, 1e-12);

    // Nearest
    amp.set_interpolation(bc::Interpolation::Nearest);
    EXPECT_NEAR(amp.evaluate(0.49), 0.0, 1e-12);
    EXPECT_NEAR(amp.evaluate(0.51), 2.0, 1e-12);
}
