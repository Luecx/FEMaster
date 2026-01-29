#include "../src/bc/support.h"
#include "../src/bc/support_collector.h"
#include "../src/cos/rectangular_system.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

using namespace fem;

// 24) Support on node region: identity vs rotated frame
TEST(BC_Support, NodeRegionIdentityAndRotated) {
    // Small model with 2 nodes
    model::Model mdl(2, 0, 0);
    mdl.set_node(0, 0.0, 0.0, 0.0);
    mdl.set_node(1, 1.0, 0.0, 0.0);

    // Identity orientation (no coordinate system)
    auto nset = std::make_shared<model::NodeRegion>("S");
    nset->add(0);
    bc::Support s_id(nset, Vec6(0, NAN, NAN, NAN, NAN, NAN));
    constraint::Equations eqs_id;
    s_id.apply(*mdl._data, eqs_id);
    ASSERT_EQ(eqs_id.size(), 1u);
    ASSERT_EQ(eqs_id[0].entries.size(), 1u);
    EXPECT_EQ(eqs_id[0].entries[0].node_id, 0);
    EXPECT_EQ(eqs_id[0].entries[0].dof, 0);
    EXPECT_NEAR(eqs_id[0].rhs, 0.0, 1e-12);

    // Rotated orientation by +90deg around Z: local x aligns with global +y
    cos::RectangularSystem rot("R", Vec3(0,1,0), Vec3(-1,0,0));
    bc::Support s_rot(nset, Vec6(0, NAN, NAN, NAN, NAN, NAN), std::make_shared<cos::RectangularSystem>(rot));
    constraint::Equations eqs_rot;
    s_rot.apply(*mdl._data, eqs_rot);
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

