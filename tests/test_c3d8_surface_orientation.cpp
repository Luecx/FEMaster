/**
 * @file test_c3d8_surface_orientation.cpp
 * @brief Verifies outward orientation of all C3D8 boundary surfaces.
 *
 * The regression specifically protects S1, whose winding was previously
 * inverted. The canonical cube makes every expected outward normal axis-aligned
 * and therefore exposes any face-ordering regression directly.
 *
 * @author Finn Eggers
 * @date 11.08.2026
 */

#include "../src/data/field.h"
#include "../src/model/solid/c3d8.h"

#include <gtest/gtest.h>

#include <array>

namespace {

using namespace fem;

TEST(C3D8SurfaceOrientation, AllFacesPointOutward) {
    const std::array<ID, 8> nodes {0, 1, 2, 3, 4, 5, 6, 7};
    model::C3D8 element(0, nodes);

    model::Field positions {
        "POSITION",
        model::FieldDomain::NODE,
        8,
        6
    };
    positions.set_zero();

    const std::array<Vec3, 8> coordinates {
        Vec3(-1, -1, -1),
        Vec3( 1, -1, -1),
        Vec3( 1,  1, -1),
        Vec3(-1,  1, -1),
        Vec3(-1, -1,  1),
        Vec3( 1, -1,  1),
        Vec3( 1,  1,  1),
        Vec3(-1,  1,  1)
    };

    for (Index node = 0; node < 8; ++node) {
        for (Dim component = 0; component < 3; ++component) {
            positions(node, component) = coordinates[static_cast<std::size_t>(node)](component);
        }
    }

    const std::array<Vec3, 6> expected {
        Vec3( 0,  0, -1), // S1
        Vec3( 0,  0,  1), // S2
        Vec3( 0, -1,  0), // S3
        Vec3( 1,  0,  0), // S4
        Vec3( 0,  1,  0), // S5
        Vec3(-1,  0,  0)  // S6
    };

    for (ID surface_id = 1; surface_id <= 6; ++surface_id) {
        const auto surface = element.surface(surface_id);
        ASSERT_NE(surface, nullptr);

        const Vec3 normal = surface->normal(positions, Vec2::Zero());
        const Vec3 reference = expected[static_cast<std::size_t>(surface_id - 1)];

        ASSERT_TRUE(normal.allFinite()) << "surface S" << surface_id;
        EXPECT_NEAR(normal(0), reference(0), 1e-12) << "surface S" << surface_id;
        EXPECT_NEAR(normal(1), reference(1), 1e-12) << "surface S" << surface_id;
        EXPECT_NEAR(normal(2), reference(2), 1e-12) << "surface S" << surface_id;
    }
}

} // namespace
