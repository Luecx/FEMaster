#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"

#include <gtest/gtest.h>

using namespace fem;

// 41) C3D8 shape functions partition of unity and Kronecker property at nodes
TEST(Elements_C3D8, ShapeFunctionsBasic) {
    model::C3D8 el(0, {0,1,2,3,4,5,6,7});
    // Sum to 1 at several points
    const Precision pts[][3] = {{0,0,0},{-0.5,0.1,0.3},{0.7,-0.3,0.2}};
    for (auto &p : pts) {
        auto N = el.shape_function(p[0], p[1], p[2]);
        Precision s = 0; for (int i=0;i<8;++i) s += N(i,0);
        EXPECT_NEAR(s, 1.0, 1e-12);
    }
    // Kronecker at corners
    const Precision corners[][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    for (int n=0;n<8;++n){
        auto N = el.shape_function(corners[n][0], corners[n][1], corners[n][2]);
        for (int i=0;i<8;++i) EXPECT_NEAR(N(i,0), i==n?1.0:0.0, 1e-12);
    }
}

// 42) C3D8 shape derivative finite-difference check
TEST(Elements_C3D8, ShapeDerivativeFD) {
    model::C3D8 el(0, {0,1,2,3,4,5,6,7});
    Precision r=0.2,s=-0.1,t=0.4, h=1e-6;
    auto dN = el.shape_derivative(r,s,t); // (8x3): dN/dr, dN/ds, dN/dt per row

    auto N0 = el.shape_function(r,s,t);
    auto Nr = el.shape_function(r+h,s,t);
    auto Ns = el.shape_function(r,s+h,t);
    auto Nt = el.shape_function(r,s,t+h);
    for (int i=0;i<8;++i){
        EXPECT_NEAR((Nr(i,0)-N0(i,0))/h, dN(i,0), 1e-5);
        EXPECT_NEAR((Ns(i,0)-N0(i,0))/h, dN(i,1), 1e-5);
        EXPECT_NEAR((Nt(i,0)-N0(i,0))/h, dN(i,2), 1e-5);
    }
}

TEST(Elements_C3D8, TopAndBottomStressAreNotFlipped) {
    model::Model model(8, 1, 0);

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = model._data->materials.activate("MAT");
    material->set_elasticity<fem::material::IsotropicElasticity>(1000.0, 0.3);
    model.solid_section("EALL", "MAT");
    model.assign_sections();

    fem::model::Field displacement("U", fem::model::FieldDomain::NODE, 8, 6);
    displacement.set_zero();
    for (int node = 0; node < 8; ++node) {
        displacement(node, 0) = (node >= 4) ? 1.0 : 0.0;
    }

    auto [stress_top, stress_bot] = model.compute_stress_top_bot(displacement, false);

    ASSERT_EQ(stress_top.rows, stress_bot.rows);
    ASSERT_EQ(stress_top.components, stress_bot.components);
    for (Index i = 0; i < stress_top.rows; ++i) {
        for (Index j = 0; j < stress_top.components; ++j) {
            EXPECT_NEAR(stress_top(i, j), stress_bot(i, j), 1e-12);
        }
    }
}
