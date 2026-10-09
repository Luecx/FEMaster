#include "../src/math/quadrature.h"

#include <gtest/gtest.h>

using namespace fem;

// 20) Quadrature weights sum (heuristic checks for canonical cases)
TEST(Math_Quadrature, WeightSums) {
    using namespace math::quadrature;
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
