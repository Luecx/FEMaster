/**
 * @file test_thermal.cpp
 * @brief Verifies scalar thermal DOF numbering and surface Robin integration.
 *
 * These focused regression tests cover the thermal infrastructure introduced by
 * the steady-state conduction path. They verify that generic DOF numbering
 * preserves arbitrary component counts, that scalar `N N^T` surface integration
 * uses the dedicated higher-order quadrature rule, and that convection produces
 * both its consistent ambient source vector and Robin boundary matrix.
 *
 * The surface tests use a planar 2 x 1 four-node rectangle so the expected
 * consistent matrices and nodal source values are available analytically.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "../src/bc/robin/convection.h"
#include "../src/mattools/numerate_dofs.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/model_data.h"

#include <gtest/gtest.h>

using namespace fem;

namespace {

// Populate a planar 2 x 1 rectangle in the global xy-plane. The same reference
// geometry is reused by the scalar surface-integration and convection tests.
void set_rectangle_positions(model::Field& positions) {
    positions(0, 0) = 0.0; positions(0, 1) = 0.0; positions(0, 2) = 0.0;
    positions(1, 0) = 2.0; positions(1, 1) = 0.0; positions(1, 2) = 0.0;
    positions(2, 0) = 2.0; positions(2, 1) = 1.0; positions(2, 2) = 0.0;
    positions(3, 0) = 0.0; positions(3, 1) = 1.0; positions(3, 2) = 0.0;
}

} // namespace

TEST(ThermalDofs, PreservesArbitraryColumnCount) {
    // Verify the scalar thermal case first: active entries must be numbered
    // contiguously while the inactive node keeps the -1 sentinel.
    SystemDofs thermal(4, 1);
    thermal << true, false, true, true;

    const auto thermal_ids = mattools::numerate_dofs(thermal);
    ASSERT_EQ(thermal_ids.rows(), 4);
    ASSERT_EQ(thermal_ids.cols(), 1);
    EXPECT_EQ(thermal_ids(0, 0), 0);
    EXPECT_EQ(thermal_ids(1, 0), -1);
    EXPECT_EQ(thermal_ids(2, 0), 1);
    EXPECT_EQ(thermal_ids(3, 0), 2);

    // Verify that generalizing the numbering routine did not alter the existing
    // six-component structural row-major numbering behavior.
    SystemDofs structural(2, 6);
    structural.fill(false);
    structural(0, 0) = true;
    structural(0, 5) = true;
    structural(1, 2) = true;

    const auto structural_ids = mattools::numerate_dofs(structural);
    ASSERT_EQ(structural_ids.rows(), 2);
    ASSERT_EQ(structural_ids.cols(), 6);
    EXPECT_EQ(structural_ids(0, 0), 0);
    EXPECT_EQ(structural_ids(0, 5), 1);
    EXPECT_EQ(structural_ids(1, 2), 2);
}

TEST(ThermalSurfaceIntegration, ShapeMatrixUsesDedicatedQuadrature) {
    // Build a bilinear rectangular surface with constant scalar coefficient one.
    // The exact integral is the standard consistent four-node surface matrix.
    model::Surface4 surface({0, 1, 2, 3});
    model::Field positions("POSITIONS", model::FieldDomain::NODE, 4, 3);
    set_rectangle_positions(positions);

    const auto matrix = surface.integrate_scalar_shape_matrix(
        positions,
        [](const Vec3&) -> Precision { return Precision(1); }
    );

    // For area A = 2, integral N N^T dGamma equals A/36 times the usual
    // [4 2 1 2; ...] consistent bilinear surface pattern.
    const Precision factor = Precision(2) / Precision(36);
    Eigen::Matrix<Precision, 4, 4> expected;
    expected << 4, 2, 1, 2,
                2, 4, 2, 1,
                1, 2, 4, 2,
                2, 1, 2, 4;
    expected *= factor;

    // Check both the complete matrix and partition-of-unity area identity.
    ASSERT_EQ(matrix.rows(), 4);
    ASSERT_EQ(matrix.cols(), 4);
    EXPECT_NEAR((matrix - expected).norm(), 0.0, 1e-12);
    EXPECT_NEAR(matrix.sum(), 2.0, 1e-12);
}

TEST(ConvectionRobin, ProducesRhsAndDirectMatrixTerms) {
    // Construct the minimal compiled-model data required by the convection
    // condition: reference nodal geometry, one surface and one surface region.
    model::ModelData data;
    data.positions_reference = std::make_shared<model::Field>(
        "POSITION_REFERENCE", model::FieldDomain::NODE, 4, 3);
    set_rectangle_positions(*data.positions_reference);

    data.surfaces.push_back(std::make_shared<model::Surface4>(
        std::array<ID, 4>{0, 1, 2, 3}));

    auto region = std::make_shared<model::SurfaceRegion>("CONVECTION_SURFACE");
    region->add(0);

    // Choose h = 2 and T_inf = 10. Over area two, the total ambient source is
    // h T_inf A = 40 and symmetry gives a nodal source of ten at every node.
    bc::Convection convection;
    convection.region_              = region;
    convection.film_coefficient_    = 2.0;
    convection.ambient_temperature_ = 10.0;

    model::Field rhs("THERMAL_LOAD", model::FieldDomain::NODE, 4, 1);
    rhs.set_zero();

    // Use an identity scalar thermal DOF map so local Robin matrix indices map
    // directly to the four rows and columns of the test matrix.
    SystemDofIds dofs(4, 1);
    dofs << 0, 1, 2, 3;

    TripletList terms;
    convection.apply(data, rhs, dofs, terms, 0.0);

    // Verify the consistently distributed ambient source contribution.
    for (Index node = 0; node < 4; ++node) {
        EXPECT_NEAR(rhs(node, 0), 10.0, 1e-12);
    }

    // Convert the assembled Robin triplets into the global sparse operator.
    SparseMatrix matrix(4, 4);
    matrix.setFromTriplets(terms.begin(), terms.end());

    // Multiplying the unit shape-product matrix by h = 2 doubles its factor.
    const Precision factor = Precision(4) / Precision(36);
    Eigen::Matrix<Precision, 4, 4> expected;
    expected << 4, 2, 1, 2,
                2, 4, 2, 1,
                1, 2, 4, 2,
                2, 1, 2, 4;
    expected *= factor;

    // Check the complete Robin operator and its integrated coefficient h A = 4.
    ASSERT_EQ(matrix.rows(), 4);
    ASSERT_EQ(matrix.cols(), 4);
    EXPECT_NEAR((DynamicMatrix(matrix) - expected).norm(), 0.0, 1e-12);
    EXPECT_NEAR(matrix.sum(), 4.0, 1e-12);
}
