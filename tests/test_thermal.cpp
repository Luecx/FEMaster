#include "../src/bc/robin/convection.h"
#include "../src/mattools/numerate_dofs.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

using namespace fem;

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
    positions << 0.0, 0.0, 0.0,
                 2.0, 0.0, 0.0,
                 2.0, 1.0, 0.0,
                 0.0, 1.0, 0.0;

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

TEST(ConvectionRobin, ProducesRhsAndSymbolicRows) {
    model::ModelData data;
    data.positions_reference = std::make_shared<model::Field>(
        "POSITION_REFERENCE", model::FieldDomain::NODE, 4, 3);
    (*data.positions_reference) << 0.0, 0.0, 0.0,
                                  2.0, 0.0, 0.0,
                                  2.0, 1.0, 0.0,
                                  0.0, 1.0, 0.0;
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
    bc::RobinEquations equations;

    convection.apply(data, rhs, equations, 0.0);

    ASSERT_EQ(equations.size(), 4u);
    for (Index node = 0; node < 4; ++node) {
        EXPECT_NEAR(rhs(node, 0), 10.0, 1e-12);
        EXPECT_EQ(equations[static_cast<std::size_t>(node)].row_node_id, node);
        EXPECT_EQ(equations[static_cast<std::size_t>(node)].row_dof, 0);

        Precision row_sum = 0;
        for (const auto& entry : equations[static_cast<std::size_t>(node)].entries) {
            row_sum += entry.coeff;
        }
        EXPECT_NEAR(row_sum, 1.0, 1e-12);
    }
}

TEST(ThermalModel, AssemblesRobinRowsCentrally) {
    model::Model model;

    SystemDofIds ids(2, 1);
    ids << 0, 1;

    bc::RobinEquations equations{
        {0, 0, {{0, 0, 2.0}, {1, 0, -1.0}}},
        {1, 0, {{0, 0, -1.0}, {1, 0, 3.0}}}
    };

    const SparseMatrix matrix = model.build_thermal_matrix(ids, equations);
    ASSERT_EQ(matrix.rows(), 2);
    ASSERT_EQ(matrix.cols(), 2);
    EXPECT_NEAR(matrix.coeff(0, 0),  2.0, 1e-12);
    EXPECT_NEAR(matrix.coeff(0, 1), -1.0, 1e-12);
    EXPECT_NEAR(matrix.coeff(1, 0), -1.0, 1e-12);
    EXPECT_NEAR(matrix.coeff(1, 1),  3.0, 1e-12);
}
