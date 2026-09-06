/**
 * @file test_bc_loads.cpp
 * @brief Tests load accumulation, point masses and inertia relief on compiled models.
 */

#include "../src/bc/load_c.h"
#include "../src/bc/load_collector.h"
#include "../src/bc/load_inertial.h"
#include "../src/loadcase/tools/inertia_relief.h"
#include "../src/model/element/point.h"
#include "../src/model/model.h"
#include "../src/section/section_point_mass.h"

#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <string>

using namespace fem;

namespace {

/**
 * Adds native POINTMASS-equivalent auxiliary point elements to all compiled nodes.
 *
 * The production POINTMASS command now represents concentrated nodal inertia by
 * post-compile PointElements carrying a shared PointMassSection. These elements
 * live in ModelData::point_elements rather than the frozen dense element array.
 */
void add_test_point_masses(model::Model& mdl, Precision mass) {
    ASSERT_TRUE(mdl._data->compiled);

    const auto nodes = mdl._data->node_sets.get(SET_NODE_ALL);
    ASSERT_NE(nodes, nullptr);

    auto section_region = std::make_shared<model::ElementRegion>("__TEST_POINTMASS");
    auto section = std::make_shared<PointMassSection>(
        section_region,
        mass,
        Vec3::Zero(),
        Vec3::Zero(),
        Vec3::Zero()
    );

    for (const ID node : *nodes) {
        const ID id = -1 - static_cast<ID>(mdl._data->point_elements.size());
        section_region->add(id);

        auto point = std::make_shared<model::PointElement>(
            id,
            std::array<ID, model::PointElement::N>{node}
        );
        point->_model_data = mdl._data.get();
        point->set_section(section);
        mdl._data->point_elements.push_back(std::move(point));
    }
}

} // namespace

TEST(BC_Loads, CLoadAdditiveOverlap) {
    model::Model mdl;
    mdl.set_node(0, 0,0,0);
    mdl.set_node(1, 1,0,0);
    mdl.compile();

    bc::LoadCollector lc("L1");

    auto r0 = std::make_shared<model::NodeRegion>("R0");
    r0->add(0);

    auto L1 = std::make_shared<bc::CLoad>();
    L1->region_ = r0;
    L1->values_ = Vec6(NAN, NAN, NAN, NAN, NAN, NAN);
    L1->values_[0] = 1.5;
    lc.add(L1);

    auto L2 = std::make_shared<bc::CLoad>();
    L2->region_ = r0;
    L2->values_ = Vec6(NAN, NAN, NAN, NAN, NAN, NAN);
    L2->values_[0] = 2.0;
    lc.add(L2);

    model::Field bc{"BC", model::FieldDomain::NODE, 2, 6};
    bc.set_zero();
    lc.apply(*mdl._data, bc, 0.0);

    EXPECT_NEAR(bc(0,0), 3.5, 1e-12);
    EXPECT_NEAR(bc(1,0), 0.0, 1e-12);
}

TEST(BC_Loads, InertialLoadIncludesPointMassesWhenEnabled) {
    model::Model mdl;
    mdl.set_node(0, 0.0, 0.0, 0.0);
    mdl.compile();

    add_test_point_masses(mdl, 2.0);

    bc::InertialLoad load;
    load.region_     = std::make_shared<model::ElementRegion>("EMPTY_REGION");
    load.center_     = Vec3::Zero();
    load.center_acc_ = Vec3(1.0, 0.0, 0.0);
    load.omega_      = Vec3::Zero();
    load.alpha_      = Vec3::Zero();

    model::Field rhs{"RHS", model::FieldDomain::NODE, 1, 6};
    rhs.set_zero();

    load.consider_point_masses_ = false;
    load.apply(*mdl._data, rhs, 0.0);
    EXPECT_NEAR(rhs(0, 0), 0.0, 1e-12);

    rhs.set_zero();
    load.consider_point_masses_ = true;
    load.apply(*mdl._data, rhs, 0.0);
    EXPECT_NEAR(rhs(0, 0), -2.0, 1e-12);
    EXPECT_NEAR(rhs(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(rhs(0, 2), 0.0, 1e-12);
}

TEST(BC_Loads, InertiaReliefBalancesPointMassOnlyModel) {
    model::Model mdl;
    mdl.set_node(0, -1.0, 0.0, 0.0);
    mdl.set_node(1,  1.0, 0.0, 0.0);
    mdl.compile();

    add_test_point_masses(mdl, 1.0);

    model::Field global_load{"GLOBAL_LOAD", model::FieldDomain::NODE, 2, 6};
    global_load.set_zero();
    global_load(0, 0) = 1.0;
    global_load(1, 0) = 1.0;

    apply_inertia_relief(*mdl._data, global_load, true);

    EXPECT_NEAR(global_load(0, 0), 0.0, 1e-10);
    EXPECT_NEAR(global_load(1, 0), 0.0, 1e-10);
}
