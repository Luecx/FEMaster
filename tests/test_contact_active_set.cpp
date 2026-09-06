/**
 * @file test_contact_active_set.cpp
 * @brief Tests current unilateral contact activation with the public assembly API.
 */

#include "../src/constraints/types/contact.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/model_data.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>

namespace {

using namespace fem;

struct ContactFixture {
    static constexpr Index Nodes = 8;

    model::ModelData             model_data;
    model::SurfaceRegion::Ptr    master_region;
    model::SurfaceRegion::Ptr    slave_region;
    std::unique_ptr<constraint::Contact> contact;
    SystemDofIds                 dofs;

    explicit ContactFixture(Precision slave_z)
        : master_region(std::make_shared<model::SurfaceRegion>("MASTER")),
          slave_region (std::make_shared<model::SurfaceRegion>("SLAVE")),
          dofs(Nodes, 6) {
        model_data.positions = std::make_shared<model::Field>(
            "POSITION", model::FieldDomain::NODE, Nodes, 3);
        model_data.positions->set_zero();

        const std::array<Vec3, 4> xy{{
            Vec3(0.0, 0.0, 0.0),
            Vec3(1.0, 0.0, 0.0),
            Vec3(1.0, 1.0, 0.0),
            Vec3(0.0, 1.0, 0.0)
        }};

        for (Index i = 0; i < 4; ++i) {
            for (Index c = 0; c < 3; ++c) {
                (*model_data.positions)(i, c) = xy[static_cast<std::size_t>(i)](c);
                (*model_data.positions)(i + 4, c) = xy[static_cast<std::size_t>(i)](c);
            }
            (*model_data.positions)(i + 4, 2) = slave_z;
        }

        model_data.surfaces.resize(2);
        model_data.surfaces[0] = std::make_shared<model::Surface4>(
            std::array<ID, 4>{0, 1, 2, 3});
        model_data.surfaces[1] = std::make_shared<model::Surface4>(
            std::array<ID, 4>{4, 5, 6, 7});

        master_region->add(0);
        slave_region->add(1);

        // Both surfaces use the same positive node ordering. Flipping the master
        // normal makes them opposing surfaces. Positive slave z then gives a
        // negative normal gap and therefore active contact.
        contact = std::make_unique<constraint::Contact>(
            master_region, slave_region,
            Precision(1000), Precision(0), true);

        dofs.setConstant(-1);
        int next = 0;
        for (Index node = 0; node < Nodes; ++node)
            for (Index component = 0; component < 3; ++component)
                dofs(node, component) = next++;
    }

    void set_slave_z(Precision z) {
        for (Index node = 4; node < 8; ++node)
            (*model_data.positions)(node, 2) = z;
    }

    Precision assemble_force_norm(std::size_t* tangent_entries = nullptr) {
        model::NodeData nodal_forces{
            "CONTACT_FORCE", model::FieldDomain::NODE, Nodes, 6};
        nodal_forces.set_zero();

        TripletList triplets;
        contact->assemble(dofs, model_data, nodal_forces, triplets);

        if (tangent_entries) *tangent_entries = triplets.size();

        Precision squared_norm = Precision(0);
        for (Index node = 0; node < Nodes; ++node)
            for (Index component = 0; component < 3; ++component)
                squared_norm += nodal_forces(node, component)
                              * nodal_forces(node, component);
        return std::sqrt(squared_norm);
    }
};

} // namespace

TEST(ContactActiveSet, PenetratingParallelSurfacesAreActive) {
    ContactFixture fixture(Precision(0.01));

    std::size_t tangent_entries = 0;
    const Precision force_norm = fixture.assemble_force_norm(&tangent_entries);

    EXPECT_GT(force_norm, Precision(0));
    EXPECT_GT(tangent_entries, 0u);
}

TEST(ContactActiveSet, SeparatedParallelSurfacesAreInactive) {
    ContactFixture fixture(Precision(-0.01));

    const Precision force_norm = fixture.assemble_force_norm();
    EXPECT_NEAR(force_norm, Precision(0), Precision(1e-12));
}
