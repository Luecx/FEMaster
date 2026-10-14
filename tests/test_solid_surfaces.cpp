/**
 * @file test_solid_surfaces.cpp
 * @brief Validates solid face numbering, interpolation nodes and outward normals.
 *
 * The expected connectivity is intentionally independent of each element's
 * surface() implementation. Canonical positive-volume geometries expose
 * inverted faces, wrong midside nodes and collapsed pyramid apex faces.
 */

#include "../src/model/solid/c3d4.h"
#include "../src/model/solid/c3d5.h"
#include "../src/model/solid/c3d6.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/solid/c3d8i.h"
#include "../src/model/solid/c3d8r.h"
#include "../src/model/solid/c3d10.h"
#include "../src/model/solid/c3d15.h"
#include "../src/model/solid/c3d20.h"
#include "../src/model/solid/c3d20r.h"

#include <gtest/gtest.h>

#include <array>
#include <numeric>
#include <vector>

using namespace fem;

namespace {

template<typename Solid, std::size_t N>
void check_faces(const std::vector<Vec3>& xyz,
                 std::size_t corner_count,
                 const std::vector<std::vector<ID>>& face_nodes) {
    ASSERT_EQ(xyz.size(), N);
    ASSERT_GE(corner_count, 4U);
    ASSERT_LE(corner_count, xyz.size());

    std::array<ID, N> ids{};
    std::iota(ids.begin(), ids.end(), ID{0});
    Solid element(0, ids);

    model::Field positions("X", model::FieldDomain::NODE, N, 3);
    for (std::size_t node = 0; node < N; ++node) {
        for (Index dim = 0; dim < 3; ++dim) {
            positions(node, dim) = xyz[node](dim);
        }
    }

    Vec3 interior = Vec3::Zero();
    for (std::size_t node = 0; node < corner_count; ++node) {
        interior += xyz[node];
    }
    interior /= static_cast<Precision>(corner_count);

    for (std::size_t f = 0; f < face_nodes.size(); ++f) {
        SCOPED_TRACE("face S" + std::to_string(f + 1));
        const auto& expected = face_nodes[f];
        auto surface = element.surface(static_cast<ID>(f + 1));
        ASSERT_NE(surface, nullptr);
        ASSERT_EQ(surface->n_nodes, static_cast<Index>(expected.size()));

        const std::vector<ID> actual(surface->nodes(), surface->nodes() + surface->n_nodes);
        EXPECT_EQ(actual, expected);

        const bool triangle = expected.size() == 3 || expected.size() == 6;
        const std::size_t corners = triangle ? 3 : 4;
        ASSERT_TRUE(expected.size() == corners || expected.size() == 2 * corners);

        // Quadratic surface mid-edge nodes must belong to the corresponding
        // volume-element edge and must appear in the surface's edge order.
        if (expected.size() == 2 * corners) {
            for (std::size_t edge = 0; edge < corners; ++edge) {
                const Vec3 midpoint =
                    (xyz[expected[edge]] + xyz[expected[(edge + 1) % corners]]) / Precision(2);
                EXPECT_NEAR((xyz[expected[corners + edge]] - midpoint).norm(),
                            Precision(0), Precision(1e-12))
                    << "mid-edge node " << (corners + edge + 1);
            }
        }

        const Vec2 local = triangle ? Vec2(Precision(1.0 / 3.0), Precision(1.0 / 3.0))
                                    : Vec2::Zero();
        const Vec3 normal = surface->normal(positions, local);
        ASSERT_TRUE(normal.allFinite());
        EXPECT_NEAR(normal.norm(), Precision(1), Precision(1e-12));
        const Vec3 on_face = surface->local_to_global(local, positions);
        EXPECT_GT(normal.dot(on_face - interior), Precision(1e-8))
            << "surface normal must point outside the solid";

        // For the affine reference geometry, the surface interpolation must
        // recover the exact planar triangle/parallelogram area.
        const Vec3 a = xyz[expected[0]];
        const Vec3 b = xyz[expected[1]];
        const Vec3 c = xyz[expected[triangle ? 2 : 3]];
        const Precision expected_area = (b - a).cross(c - a).norm()
            * (triangle ? Precision(0.5) : Precision(1));
        EXPECT_NEAR(surface->area(positions), expected_area, Precision(1e-10));
    }

    EXPECT_EQ(element.surface(static_cast<ID>(face_nodes.size() + 1)), nullptr);
}

const std::vector<Vec3> hex8 = {
    Vec3(-1, -1, -1), Vec3( 1, -1, -1), Vec3( 1,  1, -1), Vec3(-1,  1, -1),
    Vec3(-1, -1,  1), Vec3( 1, -1,  1), Vec3( 1,  1,  1), Vec3(-1,  1,  1)
};

const std::vector<Vec3> hex20 = {
    Vec3(-1, -1, -1), Vec3( 1, -1, -1), Vec3( 1,  1, -1), Vec3(-1,  1, -1),
    Vec3(-1, -1,  1), Vec3( 1, -1,  1), Vec3( 1,  1,  1), Vec3(-1,  1,  1),
    Vec3( 0, -1, -1), Vec3( 1,  0, -1), Vec3( 0,  1, -1), Vec3(-1,  0, -1),
    Vec3( 0, -1,  1), Vec3( 1,  0,  1), Vec3( 0,  1,  1), Vec3(-1,  0,  1),
    Vec3(-1, -1,  0), Vec3( 1, -1,  0), Vec3( 1,  1,  0), Vec3(-1,  1,  0)
};

const std::vector<std::vector<ID>> hex8_faces = {
    {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
    {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
};

const std::vector<std::vector<ID>> hex20_faces = {
    {0, 3, 2, 1, 11, 10, 9, 8},
    {4, 5, 6, 7, 12, 13, 14, 15},
    {0, 1, 5, 4, 8, 17, 12, 16},
    {1, 2, 6, 5, 9, 18, 13, 17},
    {2, 3, 7, 6, 10, 19, 14, 18},
    {3, 0, 4, 7, 11, 16, 15, 19}
};

} // namespace

TEST(SolidSurfaces, C3D4) {
    check_faces<model::C3D4, 4>(
        {Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,0), Vec3(0,0,1)}, 4,
        {{0,2,1}, {0,1,3}, {1,2,3}, {2,0,3}});
}

TEST(SolidSurfaces, C3D5) {
    check_faces<model::C3D5, 5>(
        {Vec3(-1,-1,-1), Vec3(1,-1,-1), Vec3(1,1,-1),
         Vec3(-1,1,-1), Vec3(0,0,1)}, 5,
        {{0,3,2,1}, {0,1,4}, {1,2,4}, {2,3,4}, {3,0,4}});
}

TEST(SolidSurfaces, C3D6) {
    check_faces<model::C3D6, 6>(
        {Vec3(1,0,-1), Vec3(0,1,-1), Vec3(0,0,-1),
         Vec3(1,0, 1), Vec3(0,1, 1), Vec3(0,0, 1)}, 6,
        {{0,2,1}, {3,4,5}, {0,1,4,3}, {1,2,5,4}, {2,0,3,5}});
}

TEST(SolidSurfaces, C3D8) {
    check_faces<model::C3D8, 8>(hex8, 8, hex8_faces);
}

TEST(SolidSurfaces, C3D8R) {
    check_faces<model::C3D8R, 8>(hex8, 8, hex8_faces);
}

TEST(SolidSurfaces, C3D8I) {
    check_faces<model::C3D8I, 8>(hex8, 8, hex8_faces);
}

TEST(SolidSurfaces, C3D10) {
    check_faces<model::C3D10, 10>(
        {Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1),
         Vec3(.5,0,0), Vec3(.5,.5,0), Vec3(0,.5,0),
         Vec3(0,0,.5), Vec3(.5,0,.5), Vec3(0,.5,.5)}, 4,
        {{0,2,1,6,5,4}, {0,1,3,4,8,7},
         {1,2,3,5,9,8}, {2,0,3,6,7,9}});
}

TEST(SolidSurfaces, C3D15) {
    check_faces<model::C3D15, 15>(
        {Vec3(0,0,-1), Vec3(1,0,-1), Vec3(0,1,-1),
         Vec3(0,0, 1), Vec3(1,0, 1), Vec3(0,1, 1),
         Vec3(.5,0,-1), Vec3(.5,.5,-1), Vec3(0,.5,-1),
         Vec3(.5,0, 1), Vec3(.5,.5, 1), Vec3(0,.5, 1),
         Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0)}, 6,
        {{0,2,1,8,7,6}, {3,4,5,9,10,11},
         {0,1,4,3,6,13,9,12}, {1,2,5,4,7,14,10,13},
         {2,0,3,5,8,12,11,14}});
}

TEST(SolidSurfaces, C3D20) {
    check_faces<model::C3D20, 20>(hex20, 8, hex20_faces);
}

TEST(SolidSurfaces, C3D20R) {
    check_faces<model::C3D20R, 20>(hex20, 8, hex20_faces);
}
