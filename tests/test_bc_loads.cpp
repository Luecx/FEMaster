/**
 * @file test_bc_loads.cpp
 * @brief Tests load accumulation, point masses and inertia relief on compiled models.
 */

#include "../src/bc/neumann/load_c.h"
#include "../src/bc/load_collector.h"
#include "../src/bc/neumann/load_inertial.h"
#include "../src/bc/neumann/load_v.h"
#include "../src/loadcase/tools/inertia_relief.h"
#include "../src/material/material.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/section/section_solid.h"

#include <gtest/gtest.h>

using namespace fem;

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

TEST(BC_Loads, VLoadDoesNotScaleByDensity) {
    model::Model mdl;
    mdl.set_node(0, 0.0, 0.0, 0.0);
    mdl.set_node(1, 1.0, 0.0, 0.0);
    mdl.set_node(2, 1.0, 1.0, 0.0);
    mdl.set_node(3, 0.0, 1.0, 0.0);
    mdl.set_node(4, 0.0, 0.0, 1.0);
    mdl.set_node(5, 1.0, 0.0, 1.0);
    mdl.set_node(6, 1.0, 1.0, 1.0);
    mdl.set_node(7, 0.0, 1.0, 1.0);
    mdl.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_density(7.0);
    mdl.add_material(material);

    const auto part = mdl._data->parts.get();
    ASSERT_NE(part, nullptr);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_   = part->elem_sets.get(SET_ELEM_ALL);
    mdl.add_section(section);
    mdl.compile();

    auto region = std::make_shared<model::ElementRegion>("SOLID");
    region->add(mdl.compiled_element_id(0));

    bc::VLoad load;
    load.region_ = region;
    load.values_ = Vec3(2.0, 0.0, 0.0);

    model::Field rhs{"RHS", model::FieldDomain::NODE, 8, 6};
    rhs.set_zero();
    load.apply(*mdl._data, rhs, 0.0);

    Precision total_x = 0.0;
    for (Index node = 0; node < 8; ++node) {
        total_x += rhs(node, 0);
    }

    // Unit volume times 2 force/volume gives a total force of 2 independent of density.
    EXPECT_NEAR(total_x, 2.0, 1e-12);
}

TEST(BC_Loads, InertialLoadScalesByDensityAndAmplitude) {
    model::Model mdl;
    mdl.set_node(0, 0.0, 0.0, 0.0);
    mdl.set_node(1, 1.0, 0.0, 0.0);
    mdl.set_node(2, 1.0, 1.0, 0.0);
    mdl.set_node(3, 0.0, 1.0, 0.0);
    mdl.set_node(4, 0.0, 0.0, 1.0);
    mdl.set_node(5, 1.0, 0.0, 1.0);
    mdl.set_node(6, 1.0, 1.0, 1.0);
    mdl.set_node(7, 0.0, 1.0, 1.0);
    mdl.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_density(7.0);
    mdl.add_material(material);

    const auto part = mdl._data->parts.get();
    ASSERT_NE(part, nullptr);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_   = part->elem_sets.get(SET_ELEM_ALL);
    mdl.add_section(section);
    mdl.compile();

    auto region = std::make_shared<model::ElementRegion>("SOLID");
    region->add(mdl.compiled_element_id(0));

    auto amplitude = std::make_shared<bc::Amplitude>("HALF", bc::Interpolation::Linear);
    amplitude->add_sample(0.0, 0.0);
    amplitude->add_sample(1.0, 0.5);

    bc::InertialLoad load;
    load.region_     = region;
    load.center_acc_ = Vec3(-2.0, 0.0, 0.0);
    load.amplitude_  = amplitude;

    model::Field rhs{"RHS", model::FieldDomain::NODE, 8, 6};
    rhs.set_zero();
    load.apply(*mdl._data, rhs, 1.0);

    Precision total_x = 0.0;
    for (Index node = 0; node < 8; ++node) {
        total_x += rhs(node, 0);
    }

    // Unit volume, density 7 and acceleration 2 give 14 force, scaled by amplitude 0.5.
    EXPECT_NEAR(total_x, 7.0, 1e-12);

    rhs.set_zero();
    load.apply(*mdl._data, rhs, 1.0, true);
    total_x = 0.0;
    for (Index node = 0; node < 8; ++node) {
        total_x += rhs(node, 0);
    }
    EXPECT_NEAR(total_x, 14.0, 1e-12);
}

TEST(BC_Loads, InertialLoadIncludesPointMassesWhenEnabled) {
    model::Model mdl;
    mdl.set_node(0, 0.0, 0.0, 0.0);
    mdl.compile();

    mdl.add_point_mass_feature(
        "NALL", 2.0, Vec3::Zero(), Vec3::Zero(), Vec3::Zero());

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

    mdl.add_point_mass_feature(
        "NALL", 1.0, Vec3::Zero(), Vec3::Zero(), Vec3::Zero());

    model::Field global_load{"GLOBAL_LOAD", model::FieldDomain::NODE, 2, 6};
    global_load.set_zero();
    global_load(0, 0) = 1.0;
    global_load(1, 0) = 1.0;

    apply_inertia_relief(*mdl._data, global_load, true);

    EXPECT_NEAR(global_load(0, 0), 0.0, 1e-10);
    EXPECT_NEAR(global_load(1, 0), 0.0, 1e-10);
}
