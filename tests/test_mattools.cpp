#include "../src/mattools/numerate_dofs.h"

#include <gtest/gtest.h>

using namespace fem;

TEST(MattoolsNumerateDofs, PreservesDynamicColumnCount) {
    SystemDofs scalar(4, 1);
    scalar << true, false, true, true;

    const auto scalar_ids = mattools::numerate_dofs(scalar);

    ASSERT_EQ(scalar_ids.rows(), 4);
    ASSERT_EQ(scalar_ids.cols(), 1);
    EXPECT_EQ(scalar_ids(0, 0), 0);
    EXPECT_EQ(scalar_ids(1, 0), -1);
    EXPECT_EQ(scalar_ids(2, 0), 1);
    EXPECT_EQ(scalar_ids(3, 0), 2);

    SystemDofs structural(2, 6);
    structural.fill(false);
    structural(0, 0) = true;
    structural(0, 5) = true;
    structural(1, 2) = true;

    const auto structural_ids = mattools::numerate_dofs(structural);

    ASSERT_EQ(structural_ids.rows(), 2);
    ASSERT_EQ(structural_ids.cols(), 6);
    EXPECT_EQ(structural_ids(0, 0), 0);
    EXPECT_EQ(structural_ids(0, 5), 1);
    EXPECT_EQ(structural_ids(1, 2), 2);
}
