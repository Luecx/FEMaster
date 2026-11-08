/**
 * @file test_model_thermal.cpp
 * @brief Tests model-level thermal system helpers.
 */

#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/truss/truss.h"

#include <gtest/gtest.h>

using namespace fem;

TEST(ModelThermalDofs, ActivatesOnlyThermalElementNodes) {
    model::Model model;

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_node(8, 2.0, 0.0, 0.0);
    model.set_node(9, 3.0, 0.0, 0.0);

    model.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);
    model.set_element<model::T3>(1, 8, 9);
    model.compile();

    const SystemDofIds ids = model.build_thermal_dof_index_matrix();

    ASSERT_EQ(ids.rows(), 10);
    ASSERT_EQ(ids.cols(), 1);

    for (Index node = 0; node < 8; ++node) {
        EXPECT_EQ(ids(node, 0), node);
    }

    EXPECT_EQ(ids(8, 0), -1);
    EXPECT_EQ(ids(9, 0), -1);
}
