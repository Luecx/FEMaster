/**
 * @file test_contact_edge_tangent.cpp
 * @brief Tests the analytic contact tangent for bounded edge projections.
 *
 * The production closest-point search is used once to create the slave-to-face
 * connectivity. The master face is then frozen while the slave surface is
 * shifted beyond one master edge. This isolates the continuous bounded edge
 * projection from contact-search topology changes and compares its analytic
 * tangent against central finite differences of the assembled residual.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#include "../src/constraints/types/contact.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/model_data.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

namespace {

using namespace fem;

struct EdgeContactFixture {
    static constexpr ID node_count = 8;

    model::ModelData data{node_count, 0, 2};
    std::shared_ptr<model::SurfaceRegion> master =
        std::make_shared<model::SurfaceRegion>("MASTER");
    std::shared_ptr<model::SurfaceRegion> slave =
        std::make_shared<model::SurfaceRegion>("SLAVE");
    SystemDofIds dofs{node_count, 6};

    EdgeContactFixture() {
        data.surfaces[0] = std::make_shared<model::Surface4>(
            std::array<ID, 4>{0, 1, 2, 3});
        data.surfaces[1] = std::make_shared<model::Surface4>(
            std::array<ID, 4>{4, 5, 6, 7});

        master->add(0);
        slave->add(1);

        data.positions = std::make_shared<model::Field>(
            "POSITION",
            model::FieldDomain::NODE,
            node_count,
            6
        );
        data.positions->set_zero();

        const std::array<Vec3, node_count> coordinates {
            Vec3(0.0, 0.0,  0.0),
            Vec3(1.0, 0.0,  0.0),
            Vec3(1.0, 1.0,  0.0),
            Vec3(0.0, 1.0,  0.0),
            Vec3(0.8, 0.4, -0.05),
            Vec3(0.9, 0.4, -0.05),
            Vec3(0.9, 0.5, -0.05),
            Vec3(0.8, 0.5, -0.05)
        };

        for (ID node = 0; node < node_count; ++node) {
            for (Dim component = 0; component < 3; ++component) {
                (*data.positions)(node, component) =
                    coordinates[static_cast<std::size_t>(node)](component);
            }
        }

        dofs.fill(-1);
        int dof = 0;
        for (ID node = 0; node < node_count; ++node) {
            for (Dim component = 0; component < 3; ++component) {
                dofs(node, component) = dof++;
            }
        }
    }

    model::NodeData forces() const {
        model::NodeData result{
            "CONTACT_FORCE",
            model::FieldDomain::NODE,
            node_count,
            6
        };
        result.set_zero();
        return result;
    }

    void shift_slave_x(Precision dx) {
        for (ID node = 4; node <= 7; ++node) {
            (*data.positions)(node, 0) += dx;
        }
    }

    DynamicVector force_vector(const model::NodeData& forces) const {
        DynamicVector result = DynamicVector::Zero(24);
        for (ID node = 0; node < node_count; ++node) {
            for (Dim component = 0; component < 3; ++component) {
                result(dofs(node, component)) = forces(node, component);
            }
        }
        return result;
    }
};

TEST(NodeSurfaceContact, BoundedEdgeTangentMatchesResidualFiniteDifference) {
    EdgeContactFixture fixture;
    constraint::Contact contact(
        fixture.master,
        fixture.slave,
        Precision(100),
        Precision(0),
        false
    );

    // Create the discrete contact connectivity while all slave nodes project
    // strictly inside the master face.
    contact.begin_update_trial();
    auto initial_forces = fixture.forces();
    TripletList initial_tangent;
    contact.assemble(fixture.dofs, fixture.data, initial_forces, initial_tangent);

    // Keep the selected master face, but move the complete slave face beyond
    // x=1. The bounded closest-point projection now lies on the master edge.
    fixture.shift_slave_x(Precision(0.2));
    contact.begin_frozen_trial();

    auto base_forces = fixture.forces();
    TripletList triplets;
    contact.assemble(fixture.dofs, fixture.data, base_forces, triplets);

    DynamicMatrix analytic = DynamicMatrix::Zero(24, 24);
    for (const auto& triplet : triplets) {
        analytic(triplet.row(), triplet.col()) += triplet.value();
    }

    constexpr Precision step = Precision(1e-7);

    // Perturb only master DOFs. The frozen contact element therefore keeps its
    // stored slave area and face while the bounded edge projection is recomputed.
    for (ID node = 0; node <= 3; ++node) {
        for (Dim component = 0; component < 3; ++component) {
            const int column = fixture.dofs(node, component);
            const Precision original = (*fixture.data.positions)(node, component);

            (*fixture.data.positions)(node, component) = original + step;
            auto plus_forces = fixture.forces();
            contact.assemble(fixture.dofs, fixture.data, plus_forces, nullptr);

            (*fixture.data.positions)(node, component) = original - step;
            auto minus_forces = fixture.forces();
            contact.assemble(fixture.dofs, fixture.data, minus_forces, nullptr);

            (*fixture.data.positions)(node, component) = original;

            const DynamicVector finite_difference =
                (fixture.force_vector(plus_forces) - fixture.force_vector(minus_forces)) /
                (Precision(2) * step);

            const Precision error = (analytic.col(column) - finite_difference).norm();
            const Precision scale = std::max(Precision(1), finite_difference.norm());

            EXPECT_LT(error / scale, Precision(2e-5))
                << "master node " << node << ", component " << component;
        }
    }

    contact.rollback_trial();
    contact.rollback_trial();
}

} // namespace
