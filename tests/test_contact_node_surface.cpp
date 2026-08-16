/**
 * @file test_contact_node_surface.cpp
 * @brief Tests frictionless node-to-surface penalty contact and its tangent.
 *
 * Parallel C3D8 faces provide analytic residual checks for separation,
 * regularized positive gap and uniform penetration. A second configuration puts
 * the slave face strictly inside the master face and compares the production
 * analytic contact tangent against finite differences of the assembled residual.
 * Two adjacent master faces additionally verify the CalculiX-style search graph,
 * frozen line-search connectivity and post-Newton topology-change detection.
 * Frozen trials also verify that retained contact elements keep both their
 * tributary slave area and their existence across the positive-gap cutoff.
 *
 * @author Finn Eggers
 * @date 12.08.2026
 */

#include "../src/constraints/types/contact.h"
#include "../src/model/geometry/surface/surface4.h"
#include "../src/model/model_data.h"
#include "../src/model/solid/c3d8.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

namespace {

using namespace fem;

constexpr Precision pi = Precision(3.141592653589793238462643383279502884L);
constexpr std::array<ID, 8> contact_nodes {4, 5, 6, 7, 8, 9, 10, 11};

struct ParallelContactFixture {
    static constexpr ID node_count = 16;

    model::ModelData data{node_count, 0, 2};

    std::shared_ptr<model::SurfaceRegion> master =
        std::make_shared<model::SurfaceRegion>("MASTER");
    std::shared_ptr<model::SurfaceRegion> slave =
        std::make_shared<model::SurfaceRegion>("SLAVE");

    SystemDofIds dofs{node_count, 6};

    explicit ParallelContactFixture(Precision slave_z) {
        const std::array<ID, 8> lower_nodes {0, 1, 2, 3, 4, 5, 6, 7};
        const std::array<ID, 8> upper_nodes {8, 9, 10, 11, 12, 13, 14, 15};

        model::C3D8 lower(0, lower_nodes);
        model::C3D8 upper(1, upper_nodes);

        data.surfaces[0] = lower.surface(2); // z = 0, outward normal +z
        data.surfaces[1] = upper.surface(1); // lower face of upper block

        master->add(0);
        slave->add(1);

        data.positions = std::make_shared<model::Field>(
            "POSITION",
            model::FieldDomain::NODE,
            node_count,
            6
        );
        data.positions->set_zero();

        const std::array<Vec3, 8> lower_coordinates {
            Vec3(0, 0, -1),
            Vec3(1, 0, -1),
            Vec3(1, 1, -1),
            Vec3(0, 1, -1),
            Vec3(0, 0,  0),
            Vec3(1, 0,  0),
            Vec3(1, 1,  0),
            Vec3(0, 1,  0)
        };

        const std::array<Vec3, 8> upper_coordinates {
            Vec3(0, 0, slave_z),
            Vec3(1, 0, slave_z),
            Vec3(1, 1, slave_z),
            Vec3(0, 1, slave_z),
            Vec3(0, 0, slave_z + 1),
            Vec3(1, 0, slave_z + 1),
            Vec3(1, 1, slave_z + 1),
            Vec3(0, 1, slave_z + 1)
        };

        for (Index local = 0; local < 8; ++local) {
            for (Dim component = 0; component < 3; ++component) {
                (*data.positions)(lower_nodes[static_cast<std::size_t>(local)], component) =
                    lower_coordinates[static_cast<std::size_t>(local)](component);
                (*data.positions)(upper_nodes[static_cast<std::size_t>(local)], component) =
                    upper_coordinates[static_cast<std::size_t>(local)](component);
            }
        }

        dofs.fill(-1);
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

    void make_slave_face_interior(Precision min_xy = Precision(0.2),
                                  Precision max_xy = Precision(0.8)) {
        const std::array<Vec2, 4> xy {
            Vec2(min_xy, min_xy),
            Vec2(max_xy, min_xy),
            Vec2(max_xy, max_xy),
            Vec2(min_xy, max_xy)
        };

        for (Index local = 0; local < 4; ++local) {
            const ID lower = 8 + local;
            const ID upper = 12 + local;
            (*data.positions)(lower, 0) = xy[static_cast<std::size_t>(local)](0);
            (*data.positions)(lower, 1) = xy[static_cast<std::size_t>(local)](1);
            (*data.positions)(upper, 0) = xy[static_cast<std::size_t>(local)](0);
            (*data.positions)(upper, 1) = xy[static_cast<std::size_t>(local)](1);
        }
    }

    void shift_slave_z(Precision dz) {
        for (ID node = 8; node <= 15; ++node) {
            (*data.positions)(node, 2) += dz;
        }
    }

    void assign_contact_dofs() {
        int dof = 0;
        for (ID node : contact_nodes) {
            for (Dim component = 0; component < 3; ++component) {
                dofs(node, component) = dof++;
            }
        }
    }
};

struct AdjacentMasterFixture {
    static constexpr ID node_count = 10;

    model::ModelData data{node_count, 0, 3};
    std::shared_ptr<model::SurfaceRegion> master = std::make_shared<model::SurfaceRegion>("MASTER");
    std::shared_ptr<model::SurfaceRegion> slave  = std::make_shared<model::SurfaceRegion>("SLAVE");
    SystemDofIds dofs{node_count, 6};

    AdjacentMasterFixture() {
        data.surfaces[0] = std::make_shared<model::Surface4>(std::array<ID, 4>{0, 1, 4, 3});
        data.surfaces[1] = std::make_shared<model::Surface4>(std::array<ID, 4>{1, 2, 5, 4});
        data.surfaces[2] = std::make_shared<model::Surface4>(std::array<ID, 4>{6, 7, 8, 9});

        master->add(0);
        master->add(1);
        slave->add(2);

        data.positions = std::make_shared<model::Field>(
            "POSITION",
            model::FieldDomain::NODE,
            node_count,
            6
        );
        data.positions->set_zero();

        const std::array<Vec3, node_count> coordinates {
            Vec3(0.0,   0.0,   0.0),
            Vec3(1.0,   0.0,   0.0),
            Vec3(2.0,   0.0,   0.0),
            Vec3(0.0,   1.0,   0.0),
            Vec3(1.0,   1.0,   0.0),
            Vec3(2.0,   1.0,   0.0),
            Vec3(1.005, 0.895, -0.1),
            Vec3(1.015, 0.895, -0.1),
            Vec3(1.015, 0.905, -0.1),
            Vec3(1.005, 0.905, -0.1)
        };

        for (ID node = 0; node < node_count; ++node) {
            for (Dim component = 0; component < 3; ++component) {
                (*data.positions)(node, component) = coordinates[static_cast<std::size_t>(node)](component);
            }
        }

        dofs.fill(-1);
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
        for (ID node = 6; node <= 9; ++node) {
            (*data.positions)(node, 0) += dx;
        }
    }
};

Precision expected_slave_force(Precision gap, Precision penalty = Precision(100)) {
    constexpr Precision slave_area = Precision(0.25);
    constexpr Precision slave_length = Precision(0.5);
    constexpr Precision smoothing_length = Precision(1e-6) * slave_length;

    const Precision factor = Precision(0.5) + std::atan(-gap / smoothing_length) / pi;
    return penalty * slave_area * gap * factor;
}

DynamicVector contact_force_vector(const ParallelContactFixture& fixture,
                                   const model::NodeData& forces) {
    DynamicVector result = DynamicVector::Zero(24);

    for (ID node : contact_nodes) {
        for (Dim component = 0; component < 3; ++component) {
            const int dof = fixture.dofs(node, component);
            if (dof >= 0) {
                result(dof) = forces(node, component);
            }
        }
    }

    return result;
}

TEST(NodeSurfaceContact, SeparatedFacesCarryNoForce) {
    ParallelContactFixture fixture(Precision(0.5));
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    auto forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, forces, nullptr);

    for (ID node = 0; node < ParallelContactFixture::node_count; ++node) {
        for (Dim component = 0; component < 6; ++component) {
            EXPECT_NEAR(forces(node, component), Precision(0), 1e-12);
        }
    }
}

TEST(NodeSurfaceContact, SmallPositiveGapUsesSmoothTensileBranch) {
    constexpr Precision gap = Precision(2.5e-7);

    ParallelContactFixture fixture(gap);
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    auto forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, forces, nullptr);

    const Precision expected = expected_slave_force(gap);
    ASSERT_GT(expected, Precision(0));

    for (ID node = 8; node <= 11; ++node) {
        EXPECT_NEAR(forces(node, 2), expected, 1e-12);
    }
}

TEST(NodeSurfaceContact, UniformPenetrationBalancesMasterAndSlaveResultants) {
    constexpr Precision gap = Precision(-0.1);

    ParallelContactFixture fixture(gap);
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    auto forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, forces, nullptr);

    const Precision expected = expected_slave_force(gap);

    Vec3 master_resultant = Vec3::Zero();
    Vec3 slave_resultant  = Vec3::Zero();

    for (ID node = 4; node <= 7; ++node) {
        master_resultant += forces.row_vec3(node);
    }
    for (ID node = 8; node <= 11; ++node) {
        slave_resultant += forces.row_vec3(node);
        EXPECT_NEAR(forces(node, 2), expected, 1e-10);
    }

    EXPECT_NEAR(slave_resultant(0), Precision(0), 1e-12);
    EXPECT_NEAR(slave_resultant(1), Precision(0), 1e-12);
    EXPECT_NEAR(slave_resultant(2), Precision(4) * expected, 1e-10);

    EXPECT_NEAR(master_resultant(0), Precision(0), 1e-12);
    EXPECT_NEAR(master_resultant(1), Precision(0), 1e-12);
    EXPECT_NEAR(master_resultant(2), -Precision(4) * expected, 1e-10);
    EXPECT_NEAR((master_resultant + slave_resultant).norm(), Precision(0), 1e-10);
}

TEST(NodeSurfaceContact, SearchGraphWalksAcrossSharedMasterEdge) {
    AdjacentMasterFixture fixture;
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    auto forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, forces, nullptr);

    // All slave nodes lie slightly to the right of x=1. Their nearest search
    // triangle centroid is still on the left face, so the CalculiX-style search
    // must cross the common x=1 edge before the complete FE-face projection.
    EXPECT_GT(forces(2, 2), Precision(0));
    EXPECT_GT(forces(5, 2), Precision(0));
    EXPECT_NEAR(forces(0, 2), Precision(0), 1e-12);
    EXPECT_NEAR(forces(3, 2), Precision(0), 1e-12);

    Vec3 total = Vec3::Zero();
    for (ID node = 0; node < AdjacentMasterFixture::node_count; ++node) {
        total += forces.row_vec3(node);
    }
    EXPECT_NEAR(total.norm(), Precision(0), 1e-10);
}

TEST(NodeSurfaceContact, FrozenTrialKeepsNewtonMasterFaceUntilCalculixRegeneration) {
    AdjacentMasterFixture fixture;
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    // Outer increment trial: the tangent-side update selects the right master
    // face and stores that discrete connectivity for all nested line-search
    // residuals.
    contact.begin_update_trial();
    auto initial_forces = fixture.forces();
    TripletList tangent;
    contact.assemble(fixture.dofs, fixture.data, initial_forces, tangent);

    EXPECT_TRUE(contact.contact_topology_changed());
    EXPECT_GT(initial_forces(2, 2), Precision(0));
    EXPECT_GT(initial_forces(5, 2), Precision(0));

    // Move the slave clearly to the left side of the shared edge. A fresh
    // search would now select face 0, but the frozen line-search transaction
    // must retain face 1 and can therefore react only through the shared edge.
    fixture.shift_slave_x(Precision(-0.03));
    contact.begin_frozen_trial();

    auto frozen_forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, frozen_forces, nullptr);

    EXPECT_NEAR(frozen_forces(0, 2), Precision(0), 1e-12);
    EXPECT_NEAR(frozen_forces(3, 2), Precision(0), 1e-12);
    EXPECT_GT(frozen_forces(1, 2), Precision(0));
    EXPECT_GT(frozen_forces(4, 2), Precision(0));

    contact.rollback_trial();

    // A later Newton regeneration is allowed to search again and moves the
    // connectivity to the left face. CalculiX iflagact depends on the NUMBER of
    // generated contact elements, so a pure face switch at unchanged count does
    // not request a discontinuity iteration.
    contact.begin_update_trial();
    auto updated_forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, updated_forces, nullptr);

    EXPECT_FALSE(contact.contact_topology_changed());
    EXPECT_GT(updated_forces(0, 2), Precision(0));
    EXPECT_GT(updated_forces(3, 2), Precision(0));

    contact.rollback_trial();
    contact.rollback_trial();
}

TEST(NodeSurfaceContact, CalculixTopologyFlagTracksContactElementCount) {
    constexpr Precision initial_gap = Precision(-1e-4);
    constexpr Precision open_gap    = Precision(1e-3);

    ParallelContactFixture fixture(initial_gap);
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    contact.begin_update_trial();

    auto initial_forces = fixture.forces();
    TripletList initial_tangent;
    contact.assemble(fixture.dofs, fixture.data, initial_forces, initial_tangent);
    EXPECT_TRUE(contact.contact_topology_changed());

    // Regenerating the same contact-element count clears iflagact.
    auto stable_forces = fixture.forces();
    TripletList stable_tangent;
    contact.assemble(fixture.dofs, fixture.data, stable_forces, stable_tangent);
    EXPECT_FALSE(contact.contact_topology_changed());

    // Moving all slave nodes beyond the generation cutoff removes the stored
    // contact elements. With the CalculiX default delcon=0.001 this count change
    // is significant and sets iflagact.
    fixture.shift_slave_z(open_gap - initial_gap);
    auto open_forces = fixture.forces();
    TripletList open_tangent;
    contact.assemble(fixture.dofs, fixture.data, open_forces, open_tangent);
    EXPECT_TRUE(contact.contact_topology_changed());

    contact.rollback_trial();
}

TEST(NodeSurfaceContact, FrozenTrialKeepsStoredSlaveArea) {
    ParallelContactFixture fixture(Precision(-0.05));
    fixture.make_slave_face_interior(Precision(0.2), Precision(0.8));
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    contact.begin_update_trial();
    auto initial_forces = fixture.forces();
    TripletList tangent;
    contact.assemble(fixture.dofs, fixture.data, initial_forces, tangent);

    const Precision initial_force = initial_forces(8, 2);
    ASSERT_LT(initial_force, Precision(0));

    // Shrink the physical slave face from area 0.36 to 0.04. The frozen trial
    // must still use the tributary area stored when the contact element was
    // generated, so its scalar slave force remains unchanged at identical gap.
    fixture.make_slave_face_interior(Precision(0.4), Precision(0.6));
    contact.begin_frozen_trial();

    auto frozen_forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, frozen_forces, nullptr);

    EXPECT_NEAR(frozen_forces(8, 2), initial_force,
                Precision(1e-10) * std::max(Precision(1), std::abs(initial_force)));

    contact.rollback_trial();
    contact.rollback_trial();
}

TEST(NodeSurfaceContact, FrozenTrialKeepsElementBeyondGenerationCutoff) {
    constexpr Precision initial_gap = Precision(-1e-4);
    constexpr Precision frozen_gap  = Precision(1e-3);

    ParallelContactFixture fixture(initial_gap);
    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    contact.begin_update_trial();
    auto initial_forces = fixture.forces();
    TripletList tangent;
    contact.assemble(fixture.dofs, fixture.data, initial_forces, tangent);
    ASSERT_LT(initial_forces(8, 2), Precision(0));

    // For A_s=0.25 the generation cutoff is 5e-4. Move the already generated
    // contact element to g=1e-3: a topology update would remove it, but a frozen
    // line-search trial must keep evaluating the same element on the smooth
    // positive-gap branch.
    fixture.shift_slave_z(frozen_gap - initial_gap);
    contact.begin_frozen_trial();

    auto frozen_forces = fixture.forces();
    contact.assemble(fixture.dofs, fixture.data, frozen_forces, nullptr);

    const Precision expected = expected_slave_force(frozen_gap);
    ASSERT_GT(expected, Precision(0));
    EXPECT_NEAR(frozen_forces(8, 2), expected, 1e-10);

    contact.rollback_trial();
    contact.rollback_trial();
}

TEST(NodeSurfaceContact, AnalyticTangentAcceptsSmallWellConditionedSurface) {
    ParallelContactFixture fixture(Precision(-0.05));
    fixture.make_slave_face_interior();
    fixture.assign_contact_dofs();

    // Uniformly scaling the complete geometry by 1e-3 reduces det(H) by 1e-12.
    // The local closest-point system remains equally well conditioned, so the
    // analytic tangent must not be rejected by an absolute determinant cutoff.
    constexpr Precision geometry_scale = Precision(1e-3);
    for (ID node = 0; node < ParallelContactFixture::node_count; ++node) {
        for (Dim component = 0; component < 3; ++component) {
            (*fixture.data.positions)(node, component) *= geometry_scale;
        }
    }

    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    auto forces = fixture.forces();
    TripletList triplets;
    contact.assemble(fixture.dofs, fixture.data, forces, triplets);

    Precision tangent_norm = Precision(0);
    for (const auto& triplet : triplets) {
        tangent_norm += std::abs(triplet.value());
    }

    EXPECT_GT(tangent_norm, Precision(0));
}

TEST(NodeSurfaceContact, AnalyticTangentMatchesResidualFiniteDifference) {
    ParallelContactFixture fixture(Precision(-0.05));
    fixture.make_slave_face_interior();
    fixture.assign_contact_dofs();

    constraint::Contact contact(fixture.master, fixture.slave, Precision(100), Precision(0), false);

    auto base_forces = fixture.forces();
    TripletList triplets;
    contact.assemble(fixture.dofs, fixture.data, base_forces, triplets);

    DynamicMatrix analytic = DynamicMatrix::Zero(24, 24);
    for (const auto& triplet : triplets) {
        analytic(triplet.row(), triplet.col()) += triplet.value();
    }

    constexpr Precision step = Precision(1e-7);

    // Master-node perturbations leave the representative slave area fixed and
    // therefore test exactly the fixed-area CalculiX contact linearization.
    for (ID node = 4; node <= 7; ++node) {
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
                (contact_force_vector(fixture, plus_forces) -
                 contact_force_vector(fixture, minus_forces)) /
                (Precision(2) * step);

            const Precision error = (analytic.col(column) - finite_difference).norm();
            const Precision scale = std::max(Precision(1), finite_difference.norm());

            EXPECT_LT(error / scale, Precision(2e-5))
                << "master node " << node << ", component " << component;
        }
    }
}

} // namespace
