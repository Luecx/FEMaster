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

TEST(Data_ModelData, ElementOffsetFieldsUpdateOutOfOrder) {
    model::Model model(4, 3, 0);

    model.set_element<model::T3>(2, 0, 1);
    ASSERT_NE(model._data->element_nodal_offsets, nullptr);
    ASSERT_NE(model._data->element_ip_offsets, nullptr);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(0), 0.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(1), 0.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(2), 0.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(3), 2.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(3), 1.0);

    model.set_element<model::T3>(0, 2, 3);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(0), 0.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(1), 2.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(2), 2.0);
    EXPECT_DOUBLE_EQ((*model._data->element_nodal_offsets)(3), 4.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(0), 0.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(1), 1.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(2), 1.0);
    EXPECT_DOUBLE_EQ((*model._data->element_ip_offsets)(3), 2.0);
}
