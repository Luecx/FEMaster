#include "../src/bc/load_collector.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

using namespace fem;

// 25) LoadCollector additive behavior on overlapping node regions (CLoad)
TEST(BC_Loads, CLoadAdditiveOverlap) {
    model::Model mdl(2, 0, 0);
    mdl.set_node(0, 0,0,0);
    mdl.set_node(1, 1,0,0);

    bc::LoadCollector lc("L1");

    // Region with node 0
    auto r0 = std::make_shared<model::NodeRegion>("R0");
    r0->add(0);

    // Two concentrated loads on same node and same dof (Fx)
    auto L1 = std::make_shared<bc::CLoad>();
    L1->region = r0;
    L1->values = Vec6(NAN, NAN, NAN, NAN, NAN, NAN);
    L1->values[0] = 1.5;
    lc.add(L1);

    auto L2 = std::make_shared<bc::CLoad>();
    L2->region = r0;
    L2->values = Vec6(NAN, NAN, NAN, NAN, NAN, NAN);
    L2->values[0] = 2.0;
    lc.add(L2);

    model::Field bc{"BC", model::FieldDomain::NODE, 2, 6};
    bc.set_zero();
    lc.apply(*mdl._data, bc, /*time=*/0.0);

    EXPECT_NEAR(bc(0,0), 3.5, 1e-12);
    EXPECT_NEAR(bc(1,0), 0.0, 1e-12);
}
