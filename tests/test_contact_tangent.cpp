/**
 * @file test_contact_tangent.cpp
 * @brief Verifies the frozen-geometry contact tangent in a rigid normal direction.
 */

#include "../src/constraints/types/contact.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/geometry/surface/surface8.h"
#include "../src/model/model_data.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <memory>

namespace {

using namespace fem;

template<typename SurfaceType>
std::array<Vec3, SurfaceType::num_nodes> planar_nodes() {
    constexpr Index N = SurfaceType::num_nodes;
    std::array<Vec3, N> nodes{};

    nodes[0] = Vec3(0.0, 0.0, 0.0);
    nodes[1] = Vec3(1.0, 0.0, 0.0);
    nodes[2] = Vec3(1.0, 1.0, 0.0);
    nodes[3] = Vec3(0.0, 1.0, 0.0);

    if constexpr (N == 8) {
        nodes[4] = Vec3(0.5, 0.0, 0.0);
        nodes[5] = Vec3(1.0, 0.5, 0.0);
        nodes[6] = Vec3(0.5, 1.0, 0.0);
        nodes[7] = Vec3(0.0, 0.5, 0.0);
    }

    return nodes;
}

template<typename SurfaceType>
class ContactTangentFixture {
public:
    static constexpr Index SurfaceNodes = SurfaceType::num_nodes;
    static constexpr Index Nodes        = 2 * SurfaceNodes;
    static constexpr Precision BaseGap  = Precision(0.01);

    model::ModelData                  model_data;
    model::SurfaceRegion::Ptr         master_region;
    model::SurfaceRegion::Ptr         slave_region;
    std::unique_ptr<constraint::Contact> contact;
    SystemDofIds                      dofs;
    Index                             n_dofs = 0;

    ContactTangentFixture()
        : master_region(std::make_shared<model::SurfaceRegion>("MASTER")),
          slave_region (std::make_shared<model::SurfaceRegion>("SLAVE")),
          dofs(Nodes, 6) {
        model_data.positions = std::make_shared<model::Field>(
            "POSITION", model::FieldDomain::NODE, Nodes, 3);
        model_data.positions->set_zero();

        const auto reference = planar_nodes<SurfaceType>();
        std::array<ID, SurfaceNodes> master_ids{};
        std::array<ID, SurfaceNodes> slave_ids{};

        for (Index i = 0; i < SurfaceNodes; ++i) {
            master_ids[static_cast<std::size_t>(i)] = i;
            slave_ids [static_cast<std::size_t>(i)] = i + SurfaceNodes;

            for (Index c = 0; c < 3; ++c) {
                (*model_data.positions)(i, c) = reference[static_cast<std::size_t>(i)](c);
                (*model_data.positions)(i + SurfaceNodes, c) =
                    reference[static_cast<std::size_t>(i)](c);
            }
            (*model_data.positions)(i + SurfaceNodes, 2) = BaseGap;
        }

        model_data.surfaces.resize(2);
        model_data.surfaces[0] = std::make_shared<SurfaceType>(master_ids);
        model_data.surfaces[1] = std::make_shared<SurfaceType>(slave_ids);

        master_region->add(0);
        slave_region->add(1);
        contact = std::make_unique<constraint::Contact>(
            master_region, slave_region,
            Precision(1000), Precision(0), true);

        dofs.setConstant(-1);
        for (Index node = 0; node < Nodes; ++node)
            for (Index component = 0; component < 3; ++component)
                dofs(node, component) = static_cast<int>(n_dofs++);
    }

    void set_slave_z(Precision z) {
        for (Index node = SurfaceNodes; node < Nodes; ++node)
            (*model_data.positions)(node, 2) = z;
    }

    DynamicVector assemble(SparseMatrix* tangent = nullptr) {
        model::NodeData nodal_forces{
            "CONTACT_FORCE", model::FieldDomain::NODE, Nodes, 6};
        nodal_forces.set_zero();

        TripletList triplets;
        contact->assemble(dofs, model_data, nodal_forces, triplets);

        DynamicVector force = DynamicVector::Zero(n_dofs);
        for (Index node = 0; node < Nodes; ++node) {
            for (Index component = 0; component < 3; ++component) {
                const int dof = dofs(node, component);
                if (dof >= 0) force(dof) = nodal_forces(node, component);
            }
        }

        if (tangent) {
            tangent->resize(n_dofs, n_dofs);
            tangent->setFromTriplets(triplets.begin(), triplets.end());
            tangent->makeCompressed();
        }

        return force;
    }

    DynamicVector rigid_slave_normal_direction() const {
        DynamicVector direction = DynamicVector::Zero(n_dofs);
        for (Index node = SurfaceNodes; node < Nodes; ++node)
            direction(dofs(node, 2)) = Precision(1);
        return direction;
    }
};

template<typename SurfaceType>
void expect_contact_tangent_matches_rigid_normal_fd() {
    ContactTangentFixture<SurfaceType> fixture;

    SparseMatrix tangent;
    fixture.set_slave_z(fixture.BaseGap);
    fixture.assemble(&tangent);

    const DynamicVector direction = fixture.rigid_slave_normal_direction();
    const DynamicVector analytic  = tangent * direction;

    const Precision h = Precision(1e-7);
    fixture.set_slave_z(fixture.BaseGap + h);
    const DynamicVector plus = fixture.assemble();
    fixture.set_slave_z(fixture.BaseGap - h);
    const DynamicVector minus = fixture.assemble();

    const DynamicVector numerical = (plus - minus) / (Precision(2) * h);
    const Precision scale = std::max({
        Precision(1), analytic.norm(), numerical.norm()
    });

    EXPECT_LT((analytic - numerical).norm() / scale, Precision(1e-5));
}

} // namespace

TEST(ContactTangent, Surface4RigidNormalDirection) {
    expect_contact_tangent_matches_rigid_normal_fd<model::Surface4>();
}

TEST(ContactTangent, Surface8RigidNormalDirection) {
    expect_contact_tangent_matches_rigid_normal_fd<model::Surface8>();
}
