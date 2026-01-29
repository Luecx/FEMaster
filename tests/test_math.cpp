#include "../src/math/quadrature.h"
#include "../src/math/interpolate.h"

#include <gtest/gtest.h>

using namespace fem;

// 20) Quadrature weights sum (heuristic checks for canonical cases)
TEST(Math_Quadrature, WeightSums) {
    using namespace quadrature;
    {
        Quadrature q(DOMAIN_ISO_LINE_A, ORDER_QUADRATIC); // [-1,1]
        Precision sumw = 0;
        for (Index i = 0; i < q.count(); ++i) sumw += q.get_point(i).w;
        EXPECT_NEAR(sumw, 2.0, 1e-12);
    }
    {
        Quadrature q(DOMAIN_ISO_QUAD, ORDER_QUADRATIC); // [-1,1]^2
        Precision sumw = 0;
        for (Index i = 0; i < q.count(); ++i) sumw += q.get_point(i).w;
        EXPECT_NEAR(sumw, 4.0, 1e-12);
    }
    {
        Quadrature q(DOMAIN_ISO_HEX, ORDER_QUADRATIC); // [-1,1]^3
        Precision sumw = 0;
        for (Index i = 0; i < q.count(); ++i) sumw += q.get_point(i).w;
        EXPECT_NEAR(sumw, 8.0, 1e-12);
    }
}

// 21) Interpolation: linear matches midpoint
TEST(Math_Interpolate, LinearMidpoint) {
    using namespace math::interpolate;

    // Create two points along x-axis with a linear field f = x
    RowMatrix xyz(2,3); xyz.setZero();
    xyz(0,0) = -1.0; xyz(1,0) = 1.0;
    RowMatrix values(2,1);
    values(0,0) = -1.0; values(1,0) = 1.0;

    Vec3 center(0.0, 0.0, 0.0);
    DynamicMatrix out = interpolate<LINEAR>(xyz, values, center);
    ASSERT_EQ(out.rows(), 1);
    ASSERT_EQ(out.cols(), 1);
    EXPECT_NEAR(out(0,0), 0.0, 1e-6);
}
