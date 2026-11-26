/**
 * @file test_constraints.cpp
 * @brief Tests constraint-system assembly, connectors and null-space mappings.
 *
 * Connector tests construct semantic nodal topology and compile it before
 * geometric constraint equations access `ModelData`. Algebraic constraint-system
 * tests remain independent of the object-model lifecycle.
 *
 * @see constraint::Connector
 * @see constraint::assemble_constraint_system
 * @see constraint::build_null_space
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 18.08.2026
 */

#include "../src/constraints/transformer/constraint_system.h"
#include "../src/constraints/transformer/constraint_transformer.h"
#include "../src/constraints/transformer/null_space.h"
#include "../src/constraints/types/connector.h"
#include "../src/constraints/types/coupling.h"
#include "../src/math/vec_util.h"
#include "../src/model/model.h"
#include "../src/cos/rectangular_system.h"

#include <gtest/gtest.h>

using namespace fem;

// Helper: small SPD matrix
static SparseMatrix eye(Index n) {
    SparseMatrix A(n,n);
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(n);
    for (Index i = 0; i < n; ++i) trips.emplace_back(i,i,1.0);
    A.setFromTriplets(trips.begin(), trips.end());
    return A;
}

// 26–28) Constraint system assembly
TEST(Constraints_Set, AssembleBasics) {
    constraint::Equations equations = {
        {{ {0, 0, 1.0} }, 0.0},                 // u0x = 0
        {{ {1, 1, 1.0}, {1, 2, 0.0} }, 5.0},    // u1y = 5
        {{ /* empty */ }, 7.0 },                // empty row
    };
    model::SystemDofs dofs(3,6); dofs.activate(0,0); dofs.activate(1,1); dofs.numerate(); // minimal map
    const auto system = constraint::assemble_constraint_system(equations, dofs, 2);
    EXPECT_EQ(system.equations, 3); // kept rows include empty unless drop-tol > 0
    EXPECT_EQ(system.dofs, 2);
    EXPECT_EQ(system.C.rows(), 3);
    EXPECT_EQ(system.C.cols(), 2);
    EXPECT_EQ(system.d.size(), 3);
}

// 27) Zero-row drop tolerance removes empty equations when enabled
TEST(Constraints_Set, ZeroRowDrop) {
    constraint::ConstraintOptions options;
    options.zero_row_drop_tolerance = 1; // enable dropping empty rows
    constraint::Equations equations = {
        {{ {0, 0, 1.0} }, 0.0}, // one real equation
        {{ /* empty */ }, 7.0 } // should be dropped
    };
    model::SystemDofs dofs(1,6); dofs.activate(0,0); dofs.numerate();
    const auto system = constraint::assemble_constraint_system(equations, dofs, 1, options);
    EXPECT_EQ(system.equations, 1);
    EXPECT_EQ(system.C.rows(), 1);
    EXPECT_EQ(system.d.size(), 1);
}

// 29) Connector equations basic
TEST(Constraints_Connector, TwoNodesSameDir) {
    // Model with positions (needed for local frames)
    model::Model mdl;
    mdl.set_node(0, 0,0,0);
    mdl.set_node(1, 1,0,0);
    mdl.compile();

    // Global x as local x
    auto cs = std::make_shared<cos::RectangularSystem>("R", Vec3(1,0,0));
    constraint::Connector conn(0, 1, cs, constraint::ConnectorType::Join); // constrain translations

    model::SystemDofs ids(2,6); for (Index n = 0; n < 2; ++n) for (Index d = 0; d < 3; ++d) ids.activate(n, d); ids.numerate();
    auto eqs = conn.get_equations(ids, *mdl._data);
    // Expect three equations (Tx,Ty,Tz) relating node 0 and 1 with opposite signs
    Index count = 0;
    for (auto& e : eqs) if (!e.entries.empty()) ++count;
    EXPECT_GE(count, 3);
}

// 30–33) Null-space rank and map transformations
TEST(Constraints_NullSpace, BuildAndTransforms) {
    // Empty set (no constraints) to exercise the identity-map path safely
    constraint::Equations equations{};
    model::SystemDofs dofs(3,6); dofs.numerate(); // 3 dofs total (n)
    const auto system = constraint::assemble_constraint_system(equations, dofs, 3);
    auto [map, rep] = constraint::build_null_space(system);
    EXPECT_TRUE(rep.feasible);
    EXPECT_EQ(map.full_size, 3);
    // Apply/recover
    DynamicVector q(map.n_master()); q.setLinSpaced(map.n_master(), 1.0, 1.0 + map.n_master() - 1);
    DynamicVector u = map.recover(q);
    EXPECT_EQ(u.size(), rep.dofs);
    DynamicVector y(rep.dofs); y.setOnes();
    DynamicVector z = map.project(y);
    EXPECT_EQ(z.size(), map.n_master());
}

// Reference node 0 with an asymmetric, non-coplanar cloud of four slave nodes
static constraint::Coupling distributing_coupling(model::Model& mdl, bool master_in_cloud = false) {
    mdl.set_node(0, 0.1, 0.2, 0.3);
    mdl.set_node(1, 1.0, 0.0, 0.0);
    mdl.set_node(2, 0.0, 2.0, 0.0);
    mdl.set_node(3, 0.0, 0.0, 1.5);
    mdl.set_node(4, 1.0, 1.0, 1.0);
    mdl.compile();

    auto slaves = std::make_shared<model::NodeRegion>("CLOUD");
    for (Index node = master_in_cloud ? 0 : 1; node <= 4; ++node) slaves->add(node);

    Dofs coupled;
    coupled.setConstant(true);
    return constraint::Coupling(0, slaves, coupled, constraint::CouplingType::STRUCTURAL);
}

static model::SystemDofs all_dofs_active(Index nodes) {
    model::SystemDofs ids(nodes, 6);
    for (Index node = 0; node < nodes; ++node)
        for (Index d = 0; d < 6; ++d) ids.activate(node, d);
    ids.numerate();
    return ids;
}

TEST(Constraints_Coupling, DistributingReproducesRigidBodyMotion) {
    model::Model mdl;
    auto         coupling = distributing_coupling(mdl);
    auto         ids      = all_dofs_active(5);
    auto         eqs      = coupling.get_equations(ids, *mdl._data);
    ASSERT_EQ(eqs.size(), 6u);

    // Rigid motion u_k = t + w x r_k of the cloud must yield q_m = [t; w]
    const Vec3 t(0.3, -0.7, 1.1);
    const Vec3 w(0.2, 0.5, -0.4);
    const auto& x = *mdl._data->positions_reference;

    DynamicMatrix q(5, 6);
    for (Index node = 0; node < 5; ++node) {
        const Vec3 u = t + w.cross(x.row_vec3(node) - x.row_vec3(0));
        q.row(node) << u(0), u(1), u(2), 0, 0, 0;
    }
    q.row(0) << t(0), t(1), t(2), w(0), w(1), w(2);

    for (const auto& eq : eqs) {
        Precision residual = 0;
        for (const auto& e : eq.entries) residual += e.coeff * q(e.node_id, e.dof);
        EXPECT_NEAR(residual, 0.0, 1e-9);
    }
}

TEST(Constraints_Coupling, DistributingTransfersMasterLoadResultant) {
    model::Model mdl;
    auto         coupling = distributing_coupling(mdl);
    auto         ids      = all_dofs_active(5);
    auto         eqs      = coupling.get_equations(ids, *mdl._data);
    ASSERT_EQ(eqs.size(), 6u);

    // Eliminating q_m = B u maps a master load b onto the slaves as f = B^T b
    const Vec6 b(3.0, 5.0, 7.0, 2.0, 4.0, 6.0);
    DynamicMatrix f = DynamicMatrix::Zero(5, 3);
    for (Index c = 0; c < 6; ++c) {
        const auto& eq = eqs[c];
        ASSERT_EQ(eq.entries.front().node_id, 0);
        ASSERT_EQ(eq.entries.front().dof,     c);
        for (Index i = 1; i < eq.entries.size(); ++i)
            f(eq.entries[i].node_id, eq.entries[i].dof) -= eq.entries[i].coeff * b(c);
    }

    const auto& x = *mdl._data->positions_reference;
    Vec3 force  = Vec3::Zero();
    Vec3 moment = Vec3::Zero();
    for (Index node = 1; node < 5; ++node) {
        const Vec3 f_k = f.row(node).transpose();
        force  += f_k;
        moment += (x.row_vec3(node) - x.row_vec3(0)).cross(f_k);
    }
    for (Index d = 0; d < 3; ++d) {
        EXPECT_NEAR(force(d),  b(d),     1e-9);
        EXPECT_NEAR(moment(d), b(d + 3), 1e-9);
    }
}

TEST(Constraints_Coupling, DistributingRejectsMasterInsideCloud) {
    model::Model mdl;
    auto         coupling = distributing_coupling(mdl, true);
    auto         ids      = all_dofs_active(5);
    EXPECT_THROW(coupling.get_equations(ids, *mdl._data), std::runtime_error);
}


// Identical support and MPC rows share the multiplier; only the support part
// is reported as reaction.
TEST(Constraints_Reactions, RedundantSupportAndMpc) {
    using Transformer = constraint::ConstraintTransformer;
    constraint::Equations equations = {
        {{{0, 0, 1.0}}, 0.0},
        {{{0, 0, 1.0}}, 0.0}
    };
    equations[0].source = constraint::EquationSourceKind::Support;
    equations[1].source = constraint::EquationSourceKind::Mpc;

    model::SystemDofs dofs(2, 6);
    dofs.activate(0, 0);
    dofs.activate(1, 0);
    dofs.numerate();

    Transformer::Options options;
    options.method = Transformer::Method::Elimination;
    Transformer transformer(equations, dofs, 2, options);
    ASSERT_EQ(transformer.unknowns(), 1);

    DynamicVector f(2);
    f << 10.0, 3.0;
    DynamicVector q(1);
    q << 3.0;

    const DynamicVector r = transformer.support_reactions(eye(2), f, q);
    EXPECT_NEAR(r[0], -5.0, 1e-6);
    EXPECT_NEAR(r[1],  0.0, 1e-6);
}


// A support recovers force along an indirect tie chain; a disconnected tie
// carries no support force.
TEST(Constraints_Reactions, TransitiveSupportComponent) {
    using Transformer = constraint::ConstraintTransformer;
    constraint::Equations equations = {
        {{{0, 0, 1.0}}, 0.0},
        {{{0, 0, -1.0}, {1, 0, 1.0}}, 0.0},
        {{{1, 0, -1.0}, {2, 0, 1.0}}, 0.0},
        {{{3, 0, -1.0}, {4, 0, 1.0}}, 0.0}
    };
    equations[0].source = constraint::EquationSourceKind::Support;
    for (Index i = 1; i < 4; ++i)
        equations[i].source = constraint::EquationSourceKind::Tie;

    model::SystemDofs dofs(5, 6);
    for (Index i = 0; i < 5; ++i) dofs.activate(i, 0);
    dofs.numerate();

    Transformer::Options options;
    options.method = Transformer::Method::Elimination;
    Transformer transformer(equations, dofs, 5, options);
    ASSERT_EQ(transformer.unknowns(), 1);

    DynamicVector f(5);
    f << 0.0, 0.0, 12.0, 5.0, -5.0;
    const DynamicVector q = DynamicVector::Zero(1);

    const DynamicVector r = transformer.support_reactions(eye(5), f, q);
    EXPECT_NEAR(r[0], -12.0, 1e-6);
    EXPECT_NEAR(r.tail(4).norm(), 0.0, 1e-6);
}
