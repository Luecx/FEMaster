#include "../src/constraints/transformer/constraint_system.h"
#include "../src/constraints/transformer/null_space.h"
#include "../src/constraints/types/connector.h"
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

// 26–28) Constraint system assembly
TEST(Constraints_Set, AssembleBasics) {
    constraint::Equations equations = {
        {{ {0, 0, 1.0} }, 0.0},                 // u0x = 0
        {{ {1, 1, 1.0}, {1, 2, 0.0} }, 5.0},    // u1y = 5
        {{ /* empty */ }, 7.0 },                // empty row
    };
    SystemDofIds dofs(3,6); dofs.setConstant(-1); dofs(0,0)=0; dofs(1,1)=1; // minimal map
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
    SystemDofIds dofs(1,6); dofs.setConstant(-1); dofs(0,0)=0;
    const auto system = constraint::assemble_constraint_system(equations, dofs, 1, options);
    EXPECT_EQ(system.equations, 1);
    EXPECT_EQ(system.C.rows(), 1);
    EXPECT_EQ(system.d.size(), 1);
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

// 30–33) Null-space rank and map transformations
TEST(Constraints_NullSpace, BuildAndTransforms) {
    // Empty set (no constraints) to exercise the identity-map path safely
    constraint::Equations equations{};
    SystemDofIds dofs(3,6); dofs.setConstant(-1); // 3 dofs total (n)
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
