/**
 * @file test_model_ids.cpp
 * @brief Tests explicit and automatic part-local surface and line identifiers.
 */

#include "../src/model/beam/b31.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"

#include <gtest/gtest.h>

#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>

using namespace fem;

namespace {

constexpr Index AUTO = model::Part::AUTO_ID_BEGIN;

void add_hex(model::Model& model, Index element_id) {
    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_element<model::C3D8>(element_id, 0, 1, 2, 3, 4, 5, 6, 7);
}

} // namespace

TEST(Model_Ids, AutomaticSurfaceIdsStartAtReservedRange) {
    model::Model model;
    add_hex(model, 0);

    model.set_surface(-1, 0, 1);
    model.set_surface(-1, 0, 2);

    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);
    EXPECT_EQ(part->surfaces.count(AUTO),     1u);
    EXPECT_EQ(part->surfaces.count(AUTO + 1), 1u);
    EXPECT_EQ(part->next_surface_id, AUTO + 2);
}

TEST(Model_Ids, MixedExplicitAndAutomaticSurfacesDoNotCollide) {
    model::Model model;
    add_hex(model, 0);

    // Automatic first, explicit afterwards, then automatic again
    model.set_surface(-1, 0, 1);
    model.set_surface( 0, 0, 2);
    model.set_surface( 5, 0, 3);
    model.set_surface(-1, 0, 4);

    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);
    EXPECT_EQ(part->surfaces.size(), 4u);
    EXPECT_EQ(part->surfaces.count(0),        1u);
    EXPECT_EQ(part->surfaces.count(5),        1u);
    EXPECT_EQ(part->surfaces.count(AUTO),     1u);
    EXPECT_EQ(part->surfaces.count(AUTO + 1), 1u);

    // Duplicate explicit ids remain an error
    EXPECT_THROW(model.set_surface(5, 0, 5), std::runtime_error);
}

TEST(Model_Ids, ExplicitIdsInReservedRangeAreRejected) {
    model::Model model;
    add_hex(model, 0);

    EXPECT_NO_THROW(model.set_surface(AUTO - 1, 0, 1));
    EXPECT_THROW   (model.set_surface(AUTO,     0, 2), std::runtime_error);
    EXPECT_THROW   (model.set_surface(AUTO + 7, 0, 3), std::runtime_error);
}

TEST(Model_Ids, AutomaticIdOverflowIsRejected) {
    model::Model model;
    add_hex(model, 0);

    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);
    part->next_surface_id = std::numeric_limits<Index>::max();

    EXPECT_THROW(model.set_surface(-1, 0, 1), std::runtime_error);
}

TEST(Model_Ids, AutomaticLineIdsUseIndependentCounter) {
    model::Model model;
    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_element<model::B31>(0, 0, 1);

    model.set_surface(-1, 0, 1);
    model.set_surface( 3, 0, 1);
    model.set_surface(-1, 0, 1);

    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);
    EXPECT_TRUE(part->surfaces.empty());
    EXPECT_EQ(part->lines.count(AUTO),     1u);
    EXPECT_EQ(part->lines.count(3),        1u);
    EXPECT_EQ(part->lines.count(AUTO + 1), 1u);
    EXPECT_EQ(part->next_line_id,    AUTO + 2);
    EXPECT_EQ(part->next_surface_id, AUTO);
}

TEST(Model_Ids, CountersAreIndependentPerPart) {
    model::Model model;

    model.add_part("A");
    add_hex(model, 0);
    model.set_surface(-1, 0, 1);
    model.set_surface(-1, 0, 2);

    model.add_part("B");
    add_hex(model, 0);
    model.set_surface(-1, 0, 1);

    const auto part_a = model._data->parts.get("A");
    const auto part_b = model._data->parts.get("B");
    ASSERT_NE(part_a, nullptr);
    ASSERT_NE(part_b, nullptr);
    EXPECT_EQ(part_a->next_surface_id, AUTO + 2);
    EXPECT_EQ(part_b->next_surface_id, AUTO + 1);
    EXPECT_EQ(part_b->surfaces.count(AUTO), 1u);
}

TEST(Model_Ids, MixedIdsCompileToDenseSurfaces) {
    model::Model model;

    model.add_part("P");
    add_hex(model, 0);
    model.set_surface(-1, 0, 1);
    model.set_surface( 7, 0, 2);
    model.set_surface(-1, 0, 3);

    model.add_instance("I1", "P", Vec3{0.0, 0.0, 0.0});
    model.add_instance("I2", "P", Vec3{2.0, 0.0, 0.0});
    model.compile();

    // Surfaces are compiled in sorted local order, explicit ids first
    ASSERT_EQ(model._data->surfaces.size(), 6u);
    EXPECT_EQ(model.compiled_surface_id(7,        "I1"), 0);
    EXPECT_EQ(model.compiled_surface_id(AUTO,     "I1"), 1);
    EXPECT_EQ(model.compiled_surface_id(AUTO + 1, "I1"), 2);
    EXPECT_EQ(model.compiled_surface_id(7,        "I2"), 3);
    EXPECT_EQ(model.compiled_surface_id(AUTO,     "I2"), 4);
    EXPECT_EQ(model.compiled_surface_id(AUTO + 1, "I2"), 5);
}

TEST(Model_Ids, LargeAutomaticSurfaceSetScalesLinearly) {
    model::Model model;

    constexpr Index n_elements = 50000;
    for (Index i = 0; i < 8; ++i) {
        model.set_node(i, static_cast<Precision>(i & 1), static_cast<Precision>((i >> 1) & 1), static_cast<Precision>(i >> 2));
    }
    for (Index i = 0; i < n_elements; ++i) {
        model.set_element<model::C3D8>(i, 0, 1, 3, 2, 4, 5, 7, 6);
    }

    const auto start = std::chrono::steady_clock::now();
    for (Index side = 1; side <= 6; ++side) {
        model.set_surface(SET_ELEM_ALL, side);
    }
    const auto elapsed = std::chrono::steady_clock::now() - start;

    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);
    EXPECT_EQ(part->surfaces.size(), 6 * n_elements);
    EXPECT_EQ(part->next_surface_id, AUTO + 6 * n_elements);

    // 300k surfaces took minutes with the quadratic search; generous bound
    EXPECT_LT(std::chrono::duration_cast<std::chrono::seconds>(elapsed).count(), 10);

    model.compile();
    EXPECT_EQ(model._data->surfaces.size(), 6 * n_elements);
}
