#include "../src/model/geometry/surface/surface4.h"   // Include path relative to the src folder
#include "../src/model/geometry/surface/surface8.h"   // Include Surface8 for testing
#include <gtest/gtest.h>

#include <array>
#include <cstdlib>

using namespace fem;

// Test fixture class for Surface4 and Surface8 using a template
template <typename SurfaceType>
class SurfaceTest : public ::testing::Test {
protected:
    SurfaceTest()
        : surface    (generate_ids()),
          node_coords("SURFACE_QUAD_TEST", model::FieldDomain::NODE, 8, 3)
    {
        // Define nodal coordinates for a square quadrilateral surface. Rows
        // 4-7 are the midside nodes used by Surface8.
        const std::array<std::array<Precision, 3>, 8> coordinates{{
            {0.0, 0.0, 0.0},
            {4.0, 0.0, 0.0},
            {4.0, 4.0, 0.0},
            {0.0, 4.0, 0.0},
            {2.0, 0.0, 0.0},
            {4.0, 2.0, 0.0},
            {2.0, 4.0, 0.0},
            {0.0, 2.0, 0.0}
        }};

        for (Index node = 0; node < 8; ++node) {
            for (Dim component = 0; component < 3; ++component) {
                node_coords(node, component) = coordinates[node][component];
            }
        }
    }

    // Helper function to generate IDs for surface nodes
    std::array<ID, SurfaceType::num_nodes> generate_ids() {
        std::array<ID, SurfaceType::num_nodes> ids{};

        for (Index i = 0; i < SurfaceType::num_nodes; ++i) {
            ids[i] = i;
        }

        return ids;
    }

    // Member variables for the fixture
    SurfaceType  surface;
    model::Field node_coords;
};

// Define the test suite name
TYPED_TEST_SUITE_P(SurfaceTest);

// Define all test cases within the test suite
TYPED_TEST_P(SurfaceTest, AreaComputation) {
    double area = this->surface.area(this->node_coords);
    double expected_area = 16.0;  // Expected area for a square surface with side length 4
    EXPECT_NEAR(area, expected_area, 1e-6);
}

TYPED_TEST_P(SurfaceTest, ProjectionWithClipping) {
    StaticVector<3> point {10.0, 10.0, 5.0};
    auto closest_point = this->surface.local_to_global(this->surface.global_to_local(point, this->node_coords, true), this->node_coords);

    EXPECT_NEAR(closest_point[0], 4.0, 1e-6);
    EXPECT_NEAR(closest_point[1], 4.0, 1e-6);
    EXPECT_NEAR(closest_point[2], 0.0, 1e-6);

    point << 2.0, 5.0, 0.0;
    closest_point = this->surface.local_to_global(this->surface.global_to_local(point, this->node_coords, true), this->node_coords);
    EXPECT_NEAR(closest_point[0], 2.0, 1e-6);
    EXPECT_NEAR(closest_point[1], 4.0, 1e-6);
    EXPECT_NEAR(closest_point[2], 0.0, 1e-6);
}

TYPED_TEST_P(SurfaceTest, ProjectionWithoutClipping) {
    StaticVector<3> point {10.0, 10.0, 5.0};
    auto closest_point = this->surface.local_to_global(this->surface.global_to_local(point, this->node_coords, false), this->node_coords);

    EXPECT_NEAR(closest_point[0], 10.0, 1e-6);
    EXPECT_NEAR(closest_point[1], 10.0, 1e-6);
    EXPECT_NEAR(closest_point[2], 0.0, 1e-6);

    point << 2.0, 5.0, 0.0;
    closest_point = this->surface.local_to_global(this->surface.global_to_local(point, this->node_coords, false), this->node_coords);
    EXPECT_NEAR(closest_point[0], 2.0, 1e-6);
    EXPECT_NEAR(closest_point[1], 5.0, 1e-6);
    EXPECT_NEAR(closest_point[2], 0.0, 1e-6);
}

TYPED_TEST_P(SurfaceTest, InBoundsSampling) {
    unsigned int seed = 12345;

    for (int i = 0; i < 10000; ++i) {
        double r = -1.0 + 2.0 * (rand_r(&seed) % 1000) / 1000.0;
        double s = -1.0 + 2.0 * (rand_r(&seed) % 1000) / 1000.0;

        bool in_bounds = this->surface.in_bounds({r, s});
        if (in_bounds) {
            EXPECT_GE(r, -1.0);
            EXPECT_LE(r, 1.0);
            EXPECT_GE(s, -1.0);
            EXPECT_LE(s, 1.0);
        }
    }
}

TYPED_TEST_P(SurfaceTest, ShapeFunctionEvaluation) {
    unsigned int seed = 23456;

    for (int i = 0; i < 100; ++i) {
        double r = -1.0 + 2.0 * (rand_r(&seed) % 1000) / 1000.0;
        double s = -1.0 + 2.0 * (rand_r(&seed) % 1000) / 1000.0;

        if (this->surface.in_bounds({r, s})) {
            auto shape_values = this->surface.shape_function(r, s);
            double sum = shape_values.sum();
            EXPECT_NEAR(sum, 1.0, 1e-6);
        }
    }
}

TYPED_TEST_P(SurfaceTest, LocalToGlobalTransformation) {
    Vec2 local_coords(-1.0, -1.0);
    auto            global_coords = this->surface.local_to_global(local_coords, this->node_coords);
    EXPECT_NEAR(global_coords[0], 0.0, 1e-6);
    EXPECT_NEAR(global_coords[1], 0.0, 1e-6);
    EXPECT_NEAR(global_coords[2], 0.0, 1e-6);

    local_coords << 1.0, -1.0;
    global_coords = this->surface.local_to_global(local_coords, this->node_coords);
    EXPECT_NEAR(global_coords[0], 4.0, 1e-6);
    EXPECT_NEAR(global_coords[1], 0.0, 1e-6);
    EXPECT_NEAR(global_coords[2], 0.0, 1e-6);

    // sample 1000 random points and make sure the mapping is fine.
    unsigned int seed = 34567;
    for (int i = 0; i < 1000; i++) {
        double r = -2.0 + 4.0 * (rand_r(&seed) % 1000) / 1000.0;
        double s = -2.0 + 4.0 * (rand_r(&seed) % 1000) / 1000.0;

        global_coords = this->surface.local_to_global({r,s}, this->node_coords);
        EXPECT_NEAR(global_coords[0], (r+1) * 2.0, 1e-6);
        EXPECT_NEAR(global_coords[1], (s+1) * 2.0, 1e-6);
        EXPECT_NEAR(global_coords[2], 0.0, 1e-6);
    }
}

// Register all the parameterized tests
REGISTER_TYPED_TEST_SUITE_P(SurfaceTest, AreaComputation, ProjectionWithClipping, ProjectionWithoutClipping, InBoundsSampling, ShapeFunctionEvaluation, LocalToGlobalTransformation);

// Instantiate the test cases for Surface4 and Surface8
typedef ::testing::Types<fem::model::Surface4, fem::model::Surface8> SurfaceTypes;
INSTANTIATE_TYPED_TEST_SUITE_P(MySurfaceTestSuite, SurfaceTest, SurfaceTypes);
