// Include necessary headers
#include "../src/model/geometry/line/line2a.h"
#include "../src/model/geometry/line/line2b.h"
#include "../src/model/geometry/line/line3a.h"
#include "../src/model/geometry/line/line3b.h"
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <random>

using namespace fem;

// Create a type list of the line types to test
typedef ::testing::Types<
    fem::model::Line2A,
    fem::model::Line2B,
    fem::model::Line3A,
    fem::model::Line3B
> LineTypes;

// Test fixture template
template <typename LineType>
class LineTest : public ::testing::Test {
protected:
    LineTest()
        : line(generateNodeIds())
    {
        // Get the number of nodes from the size of node_ids
        N = line.node_ids.size();

        // Initialize node coordinates based on the number of nodes
        node_coords.resize(N, 3);  // N nodes, 3 coordinates each (assuming 3D space)

        // For simplicity, we'll place the nodes along the x-axis from 0 to 1
        for (Index i = 0; i < N; ++i) {
            Precision x = static_cast<Precision>(i) / (N - 1);
            node_coords.row(i) << x, 0.0, 0.0;
        }
    }

    // Helper function to generate node IDs
    std::array<ID, LineType::num_nodes> generateNodeIds() {
        std::array<ID, LineType::num_nodes> ids;
        for (Index i = 0; i < LineType::num_nodes; ++i) {
            ids[i] = i;
        }
        return ids;
    }

    // Member variables for the test fixture
    LineType line;
    Eigen::Matrix<Precision, Eigen::Dynamic, 3> node_coords;
    Index N;  // Number of nodes
};

TYPED_TEST_SUITE(LineTest, LineTypes);

// A) Test that the sum of shape functions at random points is always 1
TYPED_TEST(LineTest, ShapeFunctionSumToOneAtRandomPoints) {
    Index N = this->line.node_ids.size();
    Precision min_r = this->line.min_r();
    Precision max_r = this->line.max_r();
    // Random number generator for sampling points in local coordinate range
    std::default_random_engine generator;
    std::uniform_real_distribution<Precision> distribution(min_r, max_r);

    // Test at 10 random points
    for (int i = 0; i < 10; ++i) {
        Precision r = distribution(generator);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> N_vals = this->line.shape_function(r);
        Precision sum = N_vals.sum();
        EXPECT_NEAR(sum, 1.0, 1e-6) << "Sum of shape functions at r=" << r << " is not 1.";
    }
}

// B) Test that shape functions are 1 at their respective nodes and 0 at others
TYPED_TEST(LineTest, ShapeFunctionsAreOneAtNodes) {
    Index N = this->line.node_ids.size();
    Precision min_r = this->line.min_r();
    Precision max_r = this->line.max_r();

    // Get node positions in local coordinates
    std::vector<Precision> node_positions(N);
    node_positions[0] = min_r;
    node_positions[1] = max_r;
    if (N == 3) {
        node_positions[2] = (min_r + max_r) / 2.0;
    }

    for (Index i = 0; i < N; ++i) {
        Precision r = node_positions[i];
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> N_vals = this->line.shape_function(r);
        for (Index j = 0; j < N; ++j) {
            if (i == j) {
                EXPECT_NEAR(N_vals(j), 1.0, 1e-6) << "N_" << j << "(r=" << r << ") != 1";
            } else {
                EXPECT_NEAR(N_vals(j), 0.0, 1e-6) << "N_" << j << "(r=" << r << ") != 0";
            }
        }
    }
}

// C) Test that the shape function derivatives are correct using numerical approximation
TYPED_TEST(LineTest, ShapeFunctionDerivativesAreCorrect) {
    Index N = this->line.node_ids.size();
    Precision delta = 1e-6;  // Small perturbation for numerical derivative
    Precision min_r = this->line.min_r() + delta;
    Precision max_r = this->line.max_r() - delta;
    std::default_random_engine generator;
    std::uniform_real_distribution<Precision> distribution(min_r, max_r);

    // Test at 10 random points
    for (int i = 0; i < 10; ++i) {
        Precision r = distribution(generator);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> N_plus = this->line.shape_function(r + delta);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> N_minus = this->line.shape_function(r - delta);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> numerical_derivative = (N_plus - N_minus) / (2 * delta);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> analytical_derivative = this->line.shape_derivative(r);

        for (Index j = 0; j < N; ++j) {
            EXPECT_NEAR(numerical_derivative(j), analytical_derivative(j), 1e-5)
                << "Derivative mismatch for N_" << j << " at r=" << r;
        }
    }
}

// D) Test that the second derivatives are correct using numerical approximation
TYPED_TEST(LineTest, ShapeFunctionSecondDerivativesAreCorrect) {
    Index N = this->line.node_ids.size();
    Precision delta = 1e-6;  // Small perturbation for numerical second derivative
    Precision min_r = this->line.min_r() + delta;
    Precision max_r = this->line.max_r() - delta;
    std::default_random_engine generator;
    std::uniform_real_distribution<Precision> distribution(min_r, max_r);

    // Test at 10 random points
    for (int i = 0; i < 10; ++i) {
        Precision r = distribution(generator);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> N_plus = this->line.shape_function(r + delta);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> N_ = this->line.shape_function(r);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> N_minus = this->line.shape_function(r - delta);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> numerical_second_derivative = (N_plus - 2 * N_ + N_minus) / (delta * delta);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> analytical_second_derivative = this->line.shape_second_derivative(r);

        for (Index j = 0; j < N; ++j) {
            EXPECT_NEAR(numerical_second_derivative(j), analytical_second_derivative(j), 1e-3)
                << "Second derivative mismatch for N_" << j << " at r=" << r;
        }
    }
}

// E) Test that the closest point is found correctly
TYPED_TEST(LineTest, ClosestPointIsFoundCorrectly) {
    Index N = this->line.node_ids.size();

    NodeData node_coords_system;
    node_coords_system.resize(5, 3);
    node_coords_system << 0.0, 0.0, 0.0,  // Node 0
                        1.0, 1.0, 0.0,  // Node 1
                        2.0, 0.0, 2.0,  // Node 2
                        2.5, 0.5, 0.0,  // Node 3
                        4.0, 0.0, 0.0;  // Node 4


    for(int n = 0; n < 100; n++) {
        // Define a point p in space (e.g., off the line)
        Eigen::Matrix<Precision, 3, 1> p(rand() % 1000 / 1000.0, rand() % 1000 / 1000.0, rand() % 1000 / 1000.0);


        // Find the closest point on the line in local coordinates
        Precision r_closest = this->line.global_to_local(p, node_coords_system);

        // Compute the global coordinates of the closest point
        Eigen::Matrix<Precision, 3, 1> x_closest = this->line.local_to_global(r_closest, node_coords_system);
        Precision min_distance = (x_closest - p).norm();

        // Now, check that no other point on the line is closer by subdividing the line
        Index   num_subdivisions = 1000;
        Precision min_r          = this->line.min_r();
        Precision max_r          = this->line.max_r();
        Precision delta_r        = (max_r - min_r) / (num_subdivisions - 1);

        for (Index i = 0; i < num_subdivisions; ++i) {
            Precision r = min_r + i * delta_r;
            Eigen::Matrix<Precision, 3, 1> x = this->line.local_to_global(r, node_coords_system);
            Precision distance = (x - p).norm();

            EXPECT_GE(distance, min_distance) << "Found a closer point at r=" << r;
        }
    }


}
