/**
 * @file test_data.cpp
 * @brief Tests generic fields, dictionaries and compile-time model enumeration.
 */

#include "../src/data/dict.h"
#include "../src/data/field.h"
#include "../src/data/region.h"
#include "../src/model/model.h"
#include "../src/model/truss/truss.h"

#include <cmath>
#include <gtest/gtest.h>

using namespace fem;

TEST(Data_Field, InitAndAccess) {
    model::Field field("F", model::FieldDomain::NODE, 2, 3);
    field.fill_nan();
    EXPECT_TRUE(std::isnan(field(0, 0)));

    field.set_zero();
    field(0, 0) = 1.5;
    field(0, 1) = -2.0;
    field(0, 2) = 3.0;

    auto row = field.row_vec3(0);
    EXPECT_DOUBLE_EQ(row(0), 1.5);
    EXPECT_DOUBLE_EQ(row(1), -2.0);
    EXPECT_DOUBLE_EQ(row(2), 3.0);

    Eigen::Map<RowMatrix> view(field.data(), field.rows, field.components);
    view(1, 2) = 4.25;
    EXPECT_DOUBLE_EQ(field(1, 2), 4.25);

    field.set_ones();
    EXPECT_DOUBLE_EQ(field(0, 0), 1.0);
}

TEST(Data_Dict, StringDictAndRegion) {
    model::Dict<model::NodeRegion> dict;

    EXPECT_FALSE(dict.has("A"));

    auto a = dict.create("A");
    ASSERT_TRUE(a);
    EXPECT_TRUE(dict.has("A"));
    EXPECT_EQ(a->name, std::string("A"));
    EXPECT_EQ(a->size(), 0);

    auto again = dict.activate("A");
    EXPECT_EQ(again.get(), a.get());

    dict.remove("A");
    EXPECT_FALSE(dict.has("A"));
}

TEST(Data_ModelData, SemanticDictionariesSurviveCompile) {
    model::Model model;

    model.add_part("P");
    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);

    model.set_node(10, 0.0, 0.0, 0.0);
    model.set_node(20, 1.0, 0.0, 0.0);
    model.set_element<model::T3>(30, 10, 20);
    model.add_instance("I", "P", Vec3{2.0, 0.0, 0.0});

    const auto parts_before     = model._data->parts.get("P");
    const auto instances_before = model._data->instances.get("I");

    model.compile();

    EXPECT_EQ(model._data->parts.get("P"), parts_before);
    EXPECT_EQ(model._data->instances.get("I"), instances_before);
    EXPECT_EQ(model._data->parts.get("P"), part);
    ASSERT_NE(model._data->instances.get("I"), nullptr);
    EXPECT_EQ(model._data->instances.get("I")->part, part);
    EXPECT_EQ(model.compiled_node_id(10, "I"), 0);
    EXPECT_EQ(model.compiled_node_id(20, "I"), 1);
}

TEST(Data_ModelData, ElementOffsetFieldsAreBuiltFromCompiledSparseIds) {
    model::Model model;

    model.set_node(10, 0.0, 0.0, 0.0);
    model.set_node(20, 1.0, 0.0, 0.0);
    model.set_node(30, 2.0, 0.0, 0.0);
    model.set_node(40, 3.0, 0.0, 0.0);

    model.set_element<model::T3>(20, 10, 20);
    model.set_element<model::T3>(70, 30, 40);

    EXPECT_EQ(model._data->element_nodal_offsets, nullptr);
    EXPECT_EQ(model._data->element_ip_offsets, nullptr);

    model.compile();

    ASSERT_NE(model._data->element_nodal_offsets, nullptr);
    ASSERT_NE(model._data->element_ip_offsets, nullptr);
    EXPECT_EQ(model._data->elements.size(), 2u);

    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(0), 0.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(1), 2.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(2), 4.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(0), 0.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(1), 1.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(2), 2.0);

    EXPECT_EQ(model.compiled_element_id(20), 0);
    EXPECT_EQ(model.compiled_element_id(70), 1);
}
