#include "../src/mattools/extract_connected_rows.h"
#include "../src/mattools/mask_field.h"
#include "../src/model/system_dofs.h"

#include <cmath>

#include <gtest/gtest.h>

using namespace fem;

TEST(SystemDofs, NumeratesActiveComponentsNodeMajor) {
    model::SystemDofs scalar(4, 1);
    scalar.activate(0, 0);
    scalar.activate(2, 0);
    scalar.activate(3, 0);
    scalar.numerate();

    ASSERT_EQ(scalar.nodes(), 4);
    ASSERT_EQ(scalar.width(), 1);
    ASSERT_EQ(scalar.size(), 3);
    EXPECT_EQ(scalar(0, 0), 0);
    EXPECT_EQ(scalar(1, 0), -1);
    EXPECT_EQ(scalar(2, 0), 1);
    EXPECT_EQ(scalar(3, 0), 2);

    Dofs selection;
    selection.setConstant(false);
    selection(5) = true;

    model::SystemDofs structural(2, 6);
    structural.activate(0, 0);
    structural.activate(0, selection);
    structural.activate(1, 2);
    structural.numerate();

    ASSERT_EQ(structural.nodes(), 2);
    ASSERT_EQ(structural.width(), 6);
    ASSERT_EQ(structural.size(), 3);
    EXPECT_TRUE (structural.active(0, 5));
    EXPECT_FALSE(structural.active(1, 0));
    EXPECT_EQ(structural(0, 0), 0);
    EXPECT_EQ(structural(0, 5), 1);
    EXPECT_EQ(structural(1, 2), 2);
}

TEST(SystemDofs, ReduceAndExpandRoundTrip) {
    model::SystemDofs dofs(3, 2);
    dofs.activate(0, 1);
    dofs.activate(2, 0);
    dofs.activate(2, 1);
    dofs.numerate();

    model::Field field("INPUT", model::FieldDomain::NODE, 3, 2);
    for (Index node = 0; node < 3; ++node)
        for (Index dof = 0; dof < 2; ++dof) field(node, dof) = Precision(10 * node + dof);

    const DynamicVector reduced = dofs.reduce(field);
    ASSERT_EQ(reduced.size(), 3);
    EXPECT_EQ(reduced(0), Precision(1));
    EXPECT_EQ(reduced(1), Precision(20));
    EXPECT_EQ(reduced(2), Precision(21));

    const model::Field expanded = dofs.expand(reduced);
    ASSERT_EQ(expanded.rows, 3);
    ASSERT_EQ(expanded.components, 2);
    EXPECT_EQ(expanded(0, 0), Precision(0));
    EXPECT_EQ(expanded(0, 1), Precision(1));
    EXPECT_EQ(expanded(1, 0), Precision(0));
    EXPECT_EQ(expanded(1, 1), Precision(0));
    EXPECT_EQ(expanded(2, 0), Precision(20));
    EXPECT_EQ(expanded(2, 1), Precision(21));
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


// Rows 0-1-2 form a chain over shared columns, row 3 is isolated and row 4
// only touches row 3. Seeding row 2 must select exactly the chain.
TEST(MattoolsExtractConnectedRows, SelectsTransitiveChain) {
    TripletList triplets = {
        {0, 0, 1.0},
        {1, 0, -1.0}, {1, 1, 1.0},
        {2, 1, -1.0}, {2, 2, 2.0},
        {3, 3, 1.0},
        {4, 3, -1.0}, {4, 4, 1.0}
    };
    SparseMatrix C(5, 5);
    C.setFromTriplets(triplets.begin(), triplets.end());

    std::vector<Index> rows;
    const SparseMatrix sub = mattools::extract_connected_rows(C, {2, 2}, rows);

    ASSERT_EQ(rows, (std::vector<Index>{0, 1, 2}));
    ASSERT_EQ(sub.rows(), 3);
    ASSERT_EQ(sub.cols(), 5);
    EXPECT_EQ(sub.coeff(0, 0), Precision(1));
    EXPECT_EQ(sub.coeff(1, 0), Precision(-1));
    EXPECT_EQ(sub.coeff(1, 1), Precision(1));
    EXPECT_EQ(sub.coeff(2, 1), Precision(-1));
    EXPECT_EQ(sub.coeff(2, 2), Precision(2));
    EXPECT_EQ(sub.nonZeros(), 5);

    const SparseMatrix empty = mattools::extract_connected_rows(C, {}, rows);
    EXPECT_TRUE(rows.empty());
    EXPECT_EQ(empty.rows(), 0);
}
