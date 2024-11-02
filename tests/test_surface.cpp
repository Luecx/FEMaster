#include "../src/model/geometry/surface/surface3.h"  // Include path relative to the src folder
#include <Eigen/Dense>
#include <gtest/gtest.h>

using namespace fem;

// Test fixture class for Surface3
class Surface3Test : public ::testing::Test {
    protected:
    // Use constructor initialization list to initialize Surface3 with node IDs
    Surface3Test()
        : surface({0, 1, 2})  // Initialize Surface3 with node IDs during test fixture construction
    {
        // Define node coordinates for a triangular surface
        node_coords.resize(3, 3);
        node_coords << 0.0, 0.0, 1.0,   // Node 0
            4.0, 0.0, 0.0,   // Node 1
            0.0, 1.0, 0.0;   // Node 2
    }

    // Member variables for the fixture
    fem::model::Surface3        surface;
    Eigen::Matrix<double, 3, 3> node_coords;
};

// Test for area computation
TEST_F(Surface3Test, AreaComputation) {
    double area = surface.area(node_coords);
    double expected_area = 2.8722813232690143;
    EXPECT_NEAR(area, expected_area, 1e-6);
}

// Test for projection with clipping
TEST_F(Surface3Test, ProjectionWithClipping) {
    StaticVector<3> point {3.0, 4.0, 2.0};
    auto            local_coords = surface.global_to_local(point, node_coords, true);
    // analytical result should be: 1-0.47058823529, 0.47058823529
    EXPECT_NEAR(local_coords[0], 0.5294117647, 1e-6);
    EXPECT_NEAR(local_coords[1], 0.47058823529, 1e-6);

    // pick a point which should map to a corner
    point << -2.0, 2.0, 0.0;
    local_coords = surface.global_to_local(point, node_coords, true);
    EXPECT_NEAR(local_coords[0], 0.0, 1e-6);
    EXPECT_NEAR(local_coords[1], 1.0, 1e-6);

    // pick 10 random points and make sure they are all in between 0 and 1
    for (int i = 0; i < 10; i++) {
            point << 4.0 * (rand() % 1000) / 1000.0, 1.0 * (rand() % 1000) / 1000.0, 1.0 * (rand() % 1000) / 1000.0;
            local_coords = surface.global_to_local(point, node_coords, true);
            EXPECT_GE(local_coords[0], 0.0);
            EXPECT_GE(local_coords[1], 0.0);
            EXPECT_LE(local_coords[0], 1.0);
            EXPECT_LE(local_coords[1], 1.0);
    }
}

// Test for projection without clipping
TEST_F(Surface3Test, ProjectionWithoutClipping) {
    // pick 10 random points, compute the ray from the first point to the point and check if its perpendicular to the surface
    for (int i = 0; i < 10; i++) {
        StaticVector<3> point {4.0 * (rand() % 1000) / 1000.0, 1.0 * (rand() % 1000) / 1000.0, 1.0 * (rand() % 1000) / 1000.0};
        auto            local_coords = surface.global_to_local(point, node_coords, false);
        auto            projected = surface.local_to_global(local_coords, node_coords);
        auto            diff = point - projected;
        auto diff_p1 = projected - node_coords.row(0).transpose();
        auto diff_p2 = projected - node_coords.row(1).transpose();
        auto diff_p3 = projected - node_coords.row(2).transpose();

        Precision proj_1 = diff.dot(diff_p1);
        Precision proj_2 = diff.dot(diff_p2);
        Precision proj_3 = diff.dot(diff_p3);

        EXPECT_NEAR(proj_1, 0.0, 1e-6);
        EXPECT_NEAR(proj_2, 0.0, 1e-6);
        EXPECT_NEAR(proj_3, 0.0, 1e-6);
    }
}

// Test for local to global transformation with additional points
TEST_F(Surface3Test, LocalToGlobalTransformation) {
    // Original midpoint test
    Eigen::Vector2d local_coords(0.5, 0.5);
    auto            global_coords = surface.local_to_global(local_coords, node_coords);
    EXPECT_NEAR(global_coords[0], 2.0, 1e-6);
    EXPECT_NEAR(global_coords[1], 0.5, 1e-6);
    EXPECT_NEAR(global_coords[2], 0.0, 1e-6);

    // New test cases based on the provided results
    Eigen::Vector2d local_coords1(0.5, 0.7);
    auto            global_coords1 = surface.local_to_global(local_coords1, node_coords);
    EXPECT_NEAR(global_coords1[0],  2.0, 1e-6);
    EXPECT_NEAR(global_coords1[1],  0.7, 1e-6);
    EXPECT_NEAR(global_coords1[2], -0.2, 1e-6);

    Eigen::Vector2d local_coords2(0.2, 0.4);
    auto            global_coords2 = surface.local_to_global(local_coords2, node_coords);
    StaticVector<3> expected_coords2 {2.8, 0.4, -0.1};
    EXPECT_NEAR(global_coords2[0], 0.8, 1e-6);
    EXPECT_NEAR(global_coords2[1], 0.4, 1e-6);
    EXPECT_NEAR(global_coords2[2], 0.4, 1e-6);

    Eigen::Vector2d local_coords3(0.3, 0.3);
    auto            global_coords3 = surface.local_to_global(local_coords3, node_coords);
    StaticVector<3> expected_coords3 {1.2, 0.3, 0.4};
    EXPECT_NEAR(global_coords3[0], 1.2, 1e-6);
    EXPECT_NEAR(global_coords3[1], 0.3, 1e-6);
    EXPECT_NEAR(global_coords3[2], 0.4, 1e-6);
}

// Test for shape function evaluation with more test points
TEST_F(Surface3Test, ShapeFunctionEvaluation) {
    // Test at original point (0.3, 0.3)
    Eigen::Vector2d local_coords(0.3, 0.3);
    auto            shape_values = surface.shape_function(local_coords[0], local_coords[1]);

    double sum = shape_values.sum();
    EXPECT_NEAR(sum, 1.0, 1e-6);
    EXPECT_NEAR(shape_values[0], 0.4, 1e-6);
    EXPECT_NEAR(shape_values[1], 0.3, 1e-6);
    EXPECT_NEAR(shape_values[2], 0.3, 1e-6);

    // Test at additional points r1 and r2
    Eigen::Vector2d local_coords1(0.5, 0.7);
    auto            shape_values1 = surface.shape_function(local_coords1[0], local_coords1[1]);

    sum = shape_values1.sum();
    EXPECT_NEAR(sum, 1.0, 1e-6);
    EXPECT_NEAR(shape_values1[0], -0.2, 1e-6);  // Negative values indicate outside triangle
    EXPECT_NEAR(shape_values1[1], 0.5, 1e-6);
    EXPECT_NEAR(shape_values1[2], 0.7, 1e-6);

    Eigen::Vector2d local_coords2(0.2, 0.4);
    auto            shape_values2 = surface.shape_function(local_coords2[0], local_coords2[1]);

    sum = shape_values2.sum();
    EXPECT_NEAR(sum, 1.0, 1e-6);
    EXPECT_NEAR(shape_values2[0], 0.4, 1e-6);
    EXPECT_NEAR(shape_values2[1], 0.2, 1e-6);
    EXPECT_NEAR(shape_values2[2], 0.4, 1e-6);

    // Test at triangle corners
    Eigen::Vector2d corner_coords0(0.0, 0.0);
    auto            shape_values_corner0 = surface.shape_function(corner_coords0[0], corner_coords0[1]);

    EXPECT_NEAR(shape_values_corner0[0], 1.0, 1e-6);
    EXPECT_NEAR(shape_values_corner0[1], 0.0, 1e-6);
    EXPECT_NEAR(shape_values_corner0[2], 0.0, 1e-6);

    Eigen::Vector2d corner_coords1(1.0, 0.0);
    auto            shape_values_corner1 = surface.shape_function(corner_coords1[0], corner_coords1[1]);

    EXPECT_NEAR(shape_values_corner1[0], 0.0, 1e-6);
    EXPECT_NEAR(shape_values_corner1[1], 1.0, 1e-6);
    EXPECT_NEAR(shape_values_corner1[2], 0.0, 1e-6);

    Eigen::Vector2d corner_coords2(0.0, 1.0);
    auto            shape_values_corner2 = surface.shape_function(corner_coords2[0], corner_coords2[1]);

    EXPECT_NEAR(shape_values_corner2[0], 0.0, 1e-6);
    EXPECT_NEAR(shape_values_corner2[1], 0.0, 1e-6);
    EXPECT_NEAR(shape_values_corner2[2], 1.0, 1e-6);
}
