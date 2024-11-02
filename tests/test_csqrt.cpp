//
// Created by Finn Eggers on 22.10.24.
//
#include "../src/math/csqrt.h"
#include <gtest/gtest.h>
#include "../src/core/core.h"

using namespace fem;
using namespace fem::math;

// test the csqrt function
TEST(CsqrtTest, TestCsqrt1) {
    Precision x = 3.0;
    Precision y = 4.0;
    Precision z = 5.0;
    Precision result = csqrt(x * x + y * y + z * z);
    Precision expected = 7.0710678118654755;
    EXPECT_NEAR(result, expected, 1e-10);
}

// sqrt(2) = 1.4142135623730951
TEST(CsqrtTest, TestCsqrt2) {
    Precision x = 1.0;
    Precision y = 1.0;
    Precision z = 0.0;
    Precision result = csqrt(x * x + y * y + z * z);
    Precision expected = 1.4142135623730951;
    EXPECT_NEAR(result, expected, 1e-10);
}

// sqrt(3) = 1.7320508075688772
TEST(CsqrtTest, TestCsqrt3) {
    Precision x = 1.0;
    Precision y = 1.0;
    Precision z = 1.0;
    Precision result = csqrt(x * x + y * y + z * z);
    Precision expected = 1.7320508075688772;
    EXPECT_NEAR(result, expected, 1e-10);
}

// sqrt(30 000 000) = 5477.225575051661
TEST(CsqrtTest, TestCsqrt4) {
    Precision result = csqrt(30000000);
    Precision expected = 5477.225575051661;
    EXPECT_NEAR(result, expected, 1e-10);
}