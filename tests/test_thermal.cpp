#include "../src/bc/robin/convection.h"
#include "../src/mattools/numerate_dofs.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/model_data.h"

#include <gtest/gtest.h>

using namespace fem;

namespace {

void set_rectangle_positions(model::Field& positions) {
    positions(0, 0) = 0.0; positions(0, 1) = 0.0; positions(0, 2) = 0.0;
    positions(1, 0) = 2.0; positions(1, 1) = 0.0; positions(1, 2) = 0.0;
    positions(2, 0) = 2.0; positions(2, 1) = 1.0; positions(2, 2) = 0.0;
    positions(3, 0) = 0.0; positions(3, 1) = 1.0; positions(3, 2) = 0.0;
}

} // namespace

TEST(ThermalDofs, PreservesArbitraryColumnCount) {
    SystemDofs thermal(4, 1);
    thermal << true, false, true, true;

    const auto thermal_ids = mattools::numerate_dofs(thermal);
    ASSERT_EQ(thermal_ids.rows(), 4);
    ASSERT_EQ(thermal_ids.cols(), 1);
    EXPECT_EQ(thermal_ids(0, 0), 0);
    EXPECT_EQ(thermal_ids(1, 0), -1);
    EXPECT_EQ(thermal_ids(2, 0), 1);
    EXPECT_EQ(thermal_ids(3, 0), 2);

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
    model::Surface4 surface({0, 1, 2, 3});
    model::Field positions("POSITIONS", model::FieldDomain::NODE, 4, 3);
    set_rectangle_positions(positions);

    const auto matrix = surface.integrate_scalar_shape_matrix(
        positions,
        [](const Vec3&) -> Precision { return Precision(1); }
    );

    const Precision factor = Precision(2) / Precision(36);
    Eigen::Matrix<Precision, 4, 4> expected;
    expected << 4, 2, 1, 2,
                2, 4, 2, 1,
                1, 2, 4, 2,
                2, 1, 2, 4;
    expected *= factor;

    ASSERT_EQ(matrix.rows(), 4);
    ASSERT_EQ(matrix.cols(), 4);
    EXPECT_NEAR((matrix - expected).norm(), 0.0, 1e-12);
    EXPECT_NEAR(matrix.sum(), 2.0, 1e-12);
}

TEST(ConvectionRobin, ProducesRhsAndDirectMatrixTerms) {
    model::ModelData data;
    data.positions_reference = std::make_shared<model::Field>(
        "POSITION_REFERENCE", model::FieldDomain::NODE, 4, 3);
    set_rectangle_positions(*data.positions_reference);

    data.surfaces.push_back(std::make_shared<model::Surface4>(
        std::array<ID, 4>{0, 1, 2, 3}));

    auto region = std::make_shared<model::SurfaceRegion>("CONVECTION_SURFACE");
    region->add(0);

    bc::Convection convection;
    convection.region_ = region;
    convection.film_coefficient_ = 2.0;
    convection.ambient_temperature_ = 10.0;

    model::Field rhs("THERMAL_LOAD", model::FieldDomain::NODE, 4, 1);
    rhs.set_zero();

    SystemDofIds dofs(4, 1);
    dofs << 0, 1, 2, 3;

    TripletList terms;
    convection.apply(data, rhs, dofs, terms, 0.0);

    for (Index node = 0; node < 4; ++node) {
        EXPECT_NEAR(rhs(node, 0), 10.0, 1e-12);
    }

    SparseMatrix matrix(4, 4);
    matrix.setFromTriplets(terms.begin(), terms.end());

    const Precision factor = Precision(4) / Precision(36);
    Eigen::Matrix<Precision, 4, 4> expected;
    expected << 4, 2, 1, 2,
                2, 4, 2, 1,
                1, 2, 4, 2,
                2, 1, 2, 4;
    expected *= factor;

    ASSERT_EQ(matrix.rows(), 4);
    ASSERT_EQ(matrix.cols(), 4);
    EXPECT_NEAR((DynamicMatrix(matrix) - expected).norm(), 0.0, 1e-12);
    EXPECT_NEAR(matrix.sum(), 4.0, 1e-12);
}
