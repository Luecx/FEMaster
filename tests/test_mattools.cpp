#include "../src/mattools/mask_field.h"
#include "../src/mattools/numerate_dofs.h"

#include <cmath>

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


TEST(MattoolsMaskField, KeepsSelectedEntriesAndMasksTheRest) {
    model::Field field("INPUT", model::FieldDomain::NODE, 2, 3);
    field(0, 0) = Precision(1);
    field(0, 1) = Precision(2);
    field(0, 2) = Precision(3);
    field(1, 0) = Precision(4);
    field(1, 1) = Precision(5);
    field(1, 2) = Precision(6);

    BooleanMatrix mask(2, 3);
    mask << true, false, true,
            false, true, false;

    const auto masked = mattools::mask_field(field, mask, "MASKED");

    EXPECT_EQ(masked.name, "MASKED");
    EXPECT_EQ(masked.domain, model::FieldDomain::NODE);
    EXPECT_EQ(masked(0, 0), Precision(1));
    EXPECT_TRUE(std::isnan(masked(0, 1)));
    EXPECT_EQ(masked(0, 2), Precision(3));
    EXPECT_TRUE(std::isnan(masked(1, 0)));
    EXPECT_EQ(masked(1, 1), Precision(5));
    EXPECT_TRUE(std::isnan(masked(1, 2)));
}
