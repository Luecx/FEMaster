#include "../src/mattools/numerate_dofs.h"
#include "../src/mattools/reduce_mat_to_vec.h"
#include "../src/mattools/reduce_mat_to_mat.h"
#include "../src/mattools/lump_matrix.h"
#include "../src/mattools/assemble_bc.h"

#include <gtest/gtest.h>

using namespace fem;

// 35) assemble_matrix single vs multi is heavier; covered elsewhere — here focus on numerate/reduce/lump/assemble_bc

// 36) numerate_dofs produces consecutive ids
TEST(Mattools, NumerateDofs) {
    SystemDofs mask(3,6); mask.setZero();
    mask(0,0)=true; mask(0,1)=true; mask(1,2)=true; mask(2,5)=true;
    auto ids = mattools::numerate_dofs(mask);
    // Expected order: 0,1,2,3
    std::vector<int> got;
    for (int i=0;i<ids.rows();++i) for (int j=0;j<ids.cols();++j) if (ids(i,j)>=0) got.push_back(ids(i,j));
    ASSERT_EQ(got.size(), 4u);
    EXPECT_EQ(got[0], 0);
    EXPECT_EQ(got[1], 1);
    EXPECT_EQ(got[2], 2);
    EXPECT_EQ(got[3], 3);
}

// 37) reduce_mat_to_vec picks active dofs
TEST(Mattools, ReduceMatToVec) {
    SystemDofIds ids(2,6); ids.setConstant(-1); ids(0,0)=0; ids(1,1)=1;
    model::Field loads{"LOADS", model::FieldDomain::NODE, 2, 6};
    loads.fill_nan();
    loads(0,0)=3.5; loads(1,1)=-2.0;
    DynamicVector v = mattools::reduce_mat_to_vec(ids, loads);
    ASSERT_EQ(v.size(), 2);
    EXPECT_NEAR(v(0), 3.5, 1e-12);
    EXPECT_NEAR(v(1), -2.0, 1e-12);
}

// 38) reduce_mat_to_mat submatrix extraction
TEST(Mattools, ReduceMatToMat) {
    SparseMatrix K(3,3);
    std::vector<Eigen::Triplet<Precision>> trips{{0,0,2},{0,1,1},{1,1,3},{2,2,4}};
    K.setFromTriplets(trips.begin(), trips.end());
    DynamicVector b(3); b << NAN, 0.0, NAN; // keep rows/cols 0 and 2
    SparseMatrix R = mattools::reduce_mat_to_mat(K, b);
    EXPECT_EQ(R.rows(), 2);
    EXPECT_EQ(R.cols(), 2);
}

// 39) assemble_bc ADD/SET
TEST(Mattools, AssembleBCAddSet) {
    model::Field A{"A", model::FieldDomain::NODE, 2, 6};
    model::Field B{"B", model::FieldDomain::NODE, 2, 6};
    A.fill_nan();
    B.fill_nan();
    A(0,0) = 1.0; A(1,1) = 2.0;
    B(0,0) = 3.0; B(1,1) = 2.0;
    auto C = A;
    mattools::assemble_bc(C, B, mattools::DuplicateHandling::ADD);
    EXPECT_NEAR(C(0,0), 4.0, 1e-12);
    EXPECT_NEAR(C(1,1), 4.0, 1e-12);
}

// 40) lump_matrix row sums
TEST(Mattools, LumpMatrix) {
    SparseMatrix M(3,3);
    std::vector<Eigen::Triplet<Precision>> trips{{0,0,1},{0,2,2},{1,1,3},{2,0,4}};
    M.setFromTriplets(trips.begin(), trips.end());
    SparseMatrix L = mattools::lump_matrix(M);
    EXPECT_EQ(L.rows(), 3);
    EXPECT_EQ(L.cols(), 3);
    EXPECT_NEAR(L.coeff(0,0), 3.0, 1e-12);
    EXPECT_NEAR(L.coeff(1,1), 3.0, 1e-12);
    EXPECT_NEAR(L.coeff(2,2), 4.0, 1e-12);
}
