/**
 * @file test_bc_support.cpp
 * @brief Tests nodal supports in global and rotated coordinate systems.
 *
 * Nodal topology is constructed semantically and compiled before support
 * equations consume `ModelData`, matching the production object-model lifecycle.
 *
 * @see bc::Support
 * @see model::Model::compile
 * @see cos::RectangularSystem
 *
 * @author Finn Eggers
 * @date 18.08.2026
 */

#include "../src/bc/structural/support.h"
#include "../src/bc/collector.h"
#include "../src/constraints/types/equation.h"
#include "../src/core/types_eig.h"
#include "../src/cos/rectangular_system.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

using namespace fem;

// 24) Support on node region: identity vs rotated frame
TEST(BC_Support, NodeRegionIdentityAndRotated) {
    // Small compiled model with 2 nodes
    model::Model mdl;
    mdl.set_node(0, 0.0, 0.0, 0.0);
    mdl.set_node(1, 1.0, 0.0, 0.0);
    mdl.compile();

    // Identity orientation (no coordinate system)
    auto nset = std::make_shared<model::NodeRegion>("S");
    nset->add(0);
    bc::Support s_id(nset, Vec6(0, NAN, NAN, NAN, NAN, NAN));

    model::Field          rhs{};
    constraint::Equations eqs_id{};
    SystemDofIds           system_dof_ids{};
    TripletList            matrix{};

    s_id.apply(*mdl._data, rhs, eqs_id, system_dof_ids, matrix, Precision(0), true);
    ASSERT_EQ(eqs_id.size(), 1u);
    ASSERT_EQ(eqs_id[0].entries.size(), 1u);
    EXPECT_EQ(eqs_id[0].entries[0].node_id, 0);
    EXPECT_EQ(eqs_id[0].entries[0].dof, 0);
    EXPECT_NEAR(eqs_id[0].rhs, 0.0, 1e-12);

    // Rotated orientation by +90deg around Z: local x aligns with global +y
    cos::RectangularSystem rot("R", Vec3(0,1,0), Vec3(-1,0,0));
    bc::Support s_rot(nset, Vec6(0, NAN, NAN, NAN, NAN, NAN), std::make_shared<cos::RectangularSystem>(rot));
    constraint::Equations eqs_rot{};
    s_rot.apply(*mdl._data, rhs, eqs_rot, system_dof_ids, matrix, Precision(0), true);
    // Should produce one equation with 3 entries (projection of local x onto global xyz)
    ASSERT_EQ(eqs_rot.size(), 1u);
    ASSERT_EQ(eqs_rot[0].entries.size(), 3u);
    // Projection equals [0,1,0]
    auto e0 = eqs_rot[0].entries;
    // Expect one coeff on DOF X=0 to be ~0, one on Y=1 to be ~1, one on Z=2 to be ~0
    Precision cx = 0, cy = 0, cz = 0;
    for (auto &e : e0) {
        if (e.dof == 0) cx = e.coeff;
        if (e.dof == 1) cy = e.coeff;
        if (e.dof == 2) cz = e.coeff;
        EXPECT_EQ(e.node_id, 0);
    }
    EXPECT_NEAR(cx, 0.0, 1e-12);
    EXPECT_NEAR(cy, 1.0, 1e-12);
    EXPECT_NEAR(cz, 0.0, 1e-12);
}

TEST(BC_Collectors, NamedDomainsRemainIndependent) {
    model::Model model;

    const auto loads    = model._data->load_cols.activate("SHARED");
    const auto supports = model._data->supp_cols.activate("SHARED");
    const auto thermal  = model._data->thermal_cols.activate("SHARED");

    ASSERT_NE(loads,    nullptr);
    ASSERT_NE(supports, nullptr);
    ASSERT_NE(thermal,  nullptr);

    EXPECT_NE(loads,    supports);
    EXPECT_NE(loads,    thermal);
    EXPECT_NE(supports, thermal);

    auto region = std::make_shared<model::NodeRegion>("NODES");
    auto condition = std::make_shared<bc::Support>(
        region, Vec6(0, NAN, NAN, NAN, NAN, NAN)
    );
    loads->add(condition);

    EXPECT_EQ(loads->size(),    1u);
    EXPECT_EQ(supports->size(), 0u);
    EXPECT_EQ(thermal->size(),  0u);
}
