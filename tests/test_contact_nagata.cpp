#include "../src/constraints/types/contact_nagata.h"
#include "../src/model/geometry/surface/surface3.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/geometry/surface/surface6.h"
#include "../src/model/geometry/surface/surface8.h"

#include <gtest/gtest.h>

#include <memory>
#include <vector>

using namespace fem;

namespace {

void set_position(model::Field& positions, ID node, Precision x, Precision y, Precision z = Precision(0)) {
    positions(node, 0) = x;
    positions(node, 1) = y;
    positions(node, 2) = z;
}

} // namespace

TEST(ContactNagata, TessellatesSupportedSurfacesWithoutNewVertices) {
    std::vector<model::SurfaceInterface::Ptr> surfaces(4);

    surfaces[0] = std::make_shared<model::Surface3>(std::array<ID, 3>{0, 1, 2});
    surfaces[1] = std::make_shared<model::Surface4>(std::array<ID, 4>{3, 4, 5, 6});
    surfaces[2] = std::make_shared<model::Surface6>(std::array<ID, 6>{7, 8, 9, 10, 11, 12});
    surfaces[3] = std::make_shared<model::Surface8>(std::array<ID, 8>{13, 14, 15, 16, 17, 18, 19, 20});

    model::Field positions("CONTACT_NAGATA_TEST", model::FieldDomain::NODE, 21, 3);

    set_position(positions, 0, 0, 0);
    set_position(positions, 1, 1, 0);
    set_position(positions, 2, 0, 1);

    set_position(positions, 3, 2, 0);
    set_position(positions, 4, 3, 0);
    set_position(positions, 5, 3, 1);
    set_position(positions, 6, 2, 1);

    set_position(positions, 7,  4, 0);
    set_position(positions, 8,  5, 0);
    set_position(positions, 9,  4, 1);
    set_position(positions, 10, 4.5, 0);
    set_position(positions, 11, 4.5, 0.5);
    set_position(positions, 12, 4, 0.5);

    set_position(positions, 13, 6, 0);
    set_position(positions, 14, 7, 0);
    set_position(positions, 15, 7, 1);
    set_position(positions, 16, 6, 1);
    set_position(positions, 17, 6.5, 0);
    set_position(positions, 18, 7, 0.5);
    set_position(positions, 19, 6.5, 1);
    set_position(positions, 20, 6, 0.5);

    model::SurfaceRegion region("MASTER");
    region.add(0);
    region.add(1);
    region.add(2);
    region.add(3);

    constraint::NagataContactSurface geometry(region, surfaces, positions);

    EXPECT_EQ(geometry.patch_count(), 13);

    for (const constraint::NagataPatch& patch : geometry.patches()) {
        for (ID node : patch.nodes) {
            EXPECT_GE(node, 0);
            EXPECT_LE(node, 20);
        }
    }
}

TEST(ContactNagata, PlanarTriangleReducesToLinearGeometry) {
    std::vector<model::SurfaceInterface::Ptr> surfaces(1);
    surfaces[0] = std::make_shared<model::Surface3>(std::array<ID, 3>{0, 1, 2});

    model::Field positions("CONTACT_NAGATA_TEST", model::FieldDomain::NODE, 3, 3);
    set_position(positions, 0, 0, 0);
    set_position(positions, 1, 2, 0);
    set_position(positions, 2, 0, 3);

    model::SurfaceRegion region("MASTER");
    region.add(0);

    constraint::NagataContactSurface geometry(region, surfaces, positions);
    const auto evaluation = geometry.evaluate(geometry.patch(0), Vec2(0.2, 0.3));

    EXPECT_TRUE(evaluation.valid);
    EXPECT_NEAR(evaluation.position(0), 0.4, 1e-12);
    EXPECT_NEAR(evaluation.position(1), 0.9, 1e-12);
    EXPECT_NEAR(evaluation.position(2), 0.0, 1e-12);
    EXPECT_NEAR(evaluation.normal(0), 0.0, 1e-12);
    EXPECT_NEAR(evaluation.normal(1), 0.0, 1e-12);
    EXPECT_NEAR(evaluation.normal(2), 1.0, 1e-12);

    for (const Vec3& curvature : geometry.patch(0).curvatures) {
        EXPECT_NEAR(curvature.norm(), 0.0, 1e-12);
    }
}

TEST(ContactNagata, S4PatchesUsePointerNeighbors) {
    std::vector<model::SurfaceInterface::Ptr> surfaces(1);
    surfaces[0] = std::make_shared<model::Surface4>(std::array<ID, 4>{0, 1, 2, 3});

    model::Field positions("CONTACT_NAGATA_TEST", model::FieldDomain::NODE, 4, 3);
    set_position(positions, 0, 0, 0);
    set_position(positions, 1, 1, 0);
    set_position(positions, 2, 1, 1);
    set_position(positions, 3, 0, 1);

    model::SurfaceRegion region("MASTER");
    region.add(0);

    constraint::NagataContactSurface geometry(region, surfaces, positions);

    ASSERT_EQ(geometry.patch_count(), 2);

    const constraint::NagataPatch& first  = geometry.patch(0);
    const constraint::NagataPatch& second = geometry.patch(1);

    Index shared_edges = 0;

    for (Index edge = 0; edge < 3; ++edge) {
        if (first.neighbors[static_cast<std::size_t>(edge)] == &second) {
            ++shared_edges;
            EXPECT_EQ(
                first.edge_types[static_cast<std::size_t>(edge)],
                constraint::NagataEdgeType::Internal
            );
        }
    }

    EXPECT_EQ(shared_edges, 1);
}

TEST(ContactNagata, ProjectsOntoClosedTriangle) {
    std::vector<model::SurfaceInterface::Ptr> surfaces(1);
    surfaces[0] = std::make_shared<model::Surface3>(std::array<ID, 3>{0, 1, 2});

    model::Field positions("CONTACT_NAGATA_TEST", model::FieldDomain::NODE, 3, 3);
    set_position(positions, 0, 0, 0);
    set_position(positions, 1, 1, 0);
    set_position(positions, 2, 0, 1);

    model::SurfaceRegion region("MASTER");
    region.add(0);

    constraint::NagataContactSurface geometry(region, surfaces, positions);
    const auto projection = geometry.project_on_patch(0, Vec3(0.2, 0.3, 1.0), true);

    ASSERT_TRUE(projection.valid);
    EXPECT_NEAR(projection.position(0), 0.2, 1e-9);
    EXPECT_NEAR(projection.position(1), 0.3, 1e-9);
    EXPECT_NEAR(projection.position(2), 0.0, 1e-9);
    EXPECT_NEAR(projection.distance, 1.0, 1e-9);
}
