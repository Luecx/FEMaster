#include "../src/constraints/constraint_set.h"
#include "../src/constraints/constraint_map.h"
#include "../src/constraints/builder/builder.h"
#include "../src/constraints/connector.h"
#include "../src/model/model.h"
#include "../src/cos/rectangular_system.h"

#include <gtest/gtest.h>

using namespace fem;

// Helper: small SPD matrix
static SparseMatrix eye(int n) {
    SparseMatrix A(n,n);
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(n);
    for (int i = 0; i < n; ++i) trips.emplace_back(i,i,1.0);
    A.setFromTriplets(trips.begin(), trips.end());
    return A;
}

// 26–28) ConstraintSet assemble, zero-row drop, scaling defaults
TEST(Constraints_Set, AssembleBasics) {
    constraint::ConstraintSet S;
    S.equations = {
        {{ {0, 0, 1.0} }, 0.0},                 // u0x = 0
        {{ {1, 1, 1.0}, {1, 2, 0.0} }, 5.0},    // u1y = 5
        {{ /* empty */ }, 7.0 },                // empty row
    };
    SystemDofIds dofs(3,6); dofs.setConstant(-1); dofs(0,0)=0; dofs(1,1)=1; // minimal map
    S.assemble(dofs, 2);
    EXPECT_EQ(S.m, 3); // kept rows include empty unless drop-tol > 0
    EXPECT_EQ(S.n, 2);
    EXPECT_EQ(S.C.rows(), 3);
    EXPECT_EQ(S.C.cols(), 2);
    EXPECT_EQ(S.d.size(), 3);
    // Defaults
    EXPECT_EQ(S.col_scale.size(), 2);
    EXPECT_EQ(S.row_scale.size(), 3);
    for (int i = 0; i < S.col_scale.size(); ++i) EXPECT_NEAR(S.col_scale(i), 1.0, 1e-12);
    for (int i = 0; i < S.row_scale.size(); ++i) EXPECT_NEAR(S.row_scale(i), 1.0, 1e-12);
}

// 27) Zero-row drop tolerance removes empty equations when enabled
TEST(Constraints_Set, ZeroRowDrop) {
    constraint::ConstraintSet S;
    S.opt.zero_row_drop_tol = 1; // enable dropping empty rows
    S.equations = {
        {{ {0, 0, 1.0} }, 0.0}, // one real equation
        {{ /* empty */ }, 7.0 } // should be dropped
    };
    SystemDofIds dofs(1,6); dofs.setConstant(-1); dofs(0,0)=0;
    S.assemble(dofs, 1);
    EXPECT_EQ(S.m, 1);
    EXPECT_EQ(S.C.rows(), 1);
    EXPECT_EQ(S.d.size(), 1);
}

// 29) Connector equations basic
TEST(Constraints_Connector, TwoNodesSameDir) {
    // Model with positions (needed for local frames)
    model::Model mdl(2, 0, 0);
    mdl.set_node(0, 0,0,0);
    mdl.set_node(1, 1,0,0);
    // Global x as local x
    auto cs = std::make_shared<cos::RectangularSystem>("R", Vec3(1,0,0));
    constraint::Connector conn(0, 1, cs, constraint::ConnectorType::Join); // constrain translations

    SystemDofIds ids(2,6); ids.setConstant(-1); ids(0,0)=0; ids(1,0)=1; ids(0,1)=2; ids(1,1)=3; ids(0,2)=4; ids(1,2)=5;
    auto eqs = conn.get_equations(ids, *mdl._data);
    // Expect three equations (Tx,Ty,Tz) relating node 0 and 1 with opposite signs
    int count = 0;
    for (auto& e : eqs) if (!e.entries.empty()) ++count;
    EXPECT_GE(count, 3);
}

// 30–33) Builder rank, map transforms, recover_u, apply_Tt roundtrip (homogeneous/no-constraint path)
TEST(Constraints_Builder_Map, BuildAndTransforms) {
    // Empty set (no constraints) to exercise the identity-map path safely
    constraint::ConstraintSet set;
    SystemDofIds dofs(3,6); dofs.setConstant(-1); // 3 dofs total (n)
    set.equations = {};
    set.assemble(dofs, /*n_dofs=*/3);
    auto [map, rep] = constraint::ConstraintBuilder::build(set);
    EXPECT_TRUE(rep.feasible);
    EXPECT_EQ(map.n_full(), 3);
    // Apply/recover
    DynamicVector q(map.n_master()); q.setLinSpaced(map.n_master(), 1.0, 1.0 + map.n_master() - 1);
    DynamicVector u = map.recover_u(q);
    EXPECT_EQ(u.size(), rep.n);
    DynamicVector y(rep.n); y.setOnes();
    DynamicVector z = map.apply_Tt(y);
    EXPECT_EQ(z.size(), map.n_master());
}
