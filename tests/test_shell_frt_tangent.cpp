//
// Exhaustive finite-difference verification of the MITC4 finite-rotation shell
// tangent. The test differentiates the complete nonlinear element internal-force
// vector with respect to all 24 element DOFs and compares the resulting numerical
// Jacobian against stiffness_tangent().
//

#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/shell/frt_shell_s4.h"
#include "../src/section/section_shell_integrated.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

namespace {

constexpr fem::Index NumNodes = 4;
constexpr fem::Index DofsPerNode = 6;
constexpr fem::Index NumDofs = NumNodes * DofsPerNode;

using Vec24 = fem::StaticVector<NumDofs>;
using Mat24 = fem::StaticMatrix<NumDofs, NumDofs>;

struct TangentCheckResult {
    fem::Precision h = 0.0;
    fem::Precision relative_error = std::numeric_limits<fem::Precision>::infinity();
    fem::Precision max_abs_error = std::numeric_limits<fem::Precision>::infinity();
    fem::Precision fd_asymmetry = std::numeric_limits<fem::Precision>::infinity();
    Mat24 numerical = Mat24::Zero();
};

fem::model::Model build_frt_s4_model() {
    fem::model::Model model;

    model.set_node(0,  0.00, 0.00,  0.00);
    model.set_node(1,  2.00, 0.15,  0.10);
    model.set_node(2,  2.20, 1.35,  0.32);
    model.set_node(3, -0.10, 1.10, -0.06);
    model.set_element<fem::model::FRTShellS4>(0, 0, 1, 2, 3);

    auto material = std::make_shared<fem::material::Material>("MAT");
    material->set_elasticity<fem::material::IsotropicElasticity>(210000.0, 0.30);
    model.add_material(material);

    model.add_section(std::make_shared<fem::IntegratedShellSection>(
        material,
        model._data->parts.get()->elem_sets.get(SET_ELEM_ALL),
        0.08
    ));

    model.compile();
    return model;
}

fem::model::Field displacement_from_vector(const Vec24& q) {
    fem::model::Field displacement{
        "U", fem::model::FieldDomain::NODE, NumNodes, DofsPerNode
    };
    displacement.set_zero();

    for (fem::Index node = 0; node < NumNodes; ++node)
        for (fem::Index dof = 0; dof < DofsPerNode; ++dof)
            displacement(node, dof) = q(DofsPerNode * node + dof);

    return displacement;
}

Vec24 gather_force(const fem::model::NodeData& nodal_forces) {
    Vec24 force = Vec24::Zero();
    for (fem::Index node = 0; node < NumNodes; ++node)
        for (fem::Index dof = 0; dof < DofsPerNode; ++dof)
            force(DofsPerNode * node + dof) = nodal_forces(node, dof);
    return force;
}

Vec24 internal_force(fem::model::FRTShellS4& element, const Vec24& q) {
    fem::model::NodeData nodal_forces{
        "INTERNAL_FORCES", fem::model::FieldDomain::NODE, NumNodes, DofsPerNode
    };
    nodal_forces.set_zero();

    std::array<fem::Precision, NumDofs * NumDofs> storage{};
    const fem::model::Field displacement = displacement_from_vector(q);
    element.stiffness_tangent(storage.data(), nodal_forces, displacement);

    return gather_force(nodal_forces);
}

Mat24 analytic_tangent(fem::model::FRTShellS4& element, const Vec24& q) {
    fem::model::NodeData nodal_forces{
        "INTERNAL_FORCES", fem::model::FieldDomain::NODE, NumNodes, DofsPerNode
    };
    nodal_forces.set_zero();

    std::array<fem::Precision, NumDofs * NumDofs> storage{};
    const fem::model::Field displacement = displacement_from_vector(q);
    const fem::DynamicMatrix mapped = element.stiffness_tangent(
        storage.data(), nodal_forces, displacement
    );

    Mat24 tangent = mapped;
    return tangent;
}

Mat24 finite_difference_tangent_5point(
    fem::model::FRTShellS4& element,
    const Vec24& q,
    fem::Precision h
) {
    Mat24 numerical = Mat24::Zero();

    for (fem::Index column = 0; column < NumDofs; ++column) {
        const fem::Precision step = h * std::max(fem::Precision(1), std::abs(q(column)));

        Vec24 qm2 = q;
        Vec24 qm1 = q;
        Vec24 qp1 = q;
        Vec24 qp2 = q;
        qm2(column) -= fem::Precision(2) * step;
        qm1(column) -= step;
        qp1(column) += step;
        qp2(column) += fem::Precision(2) * step;

        const Vec24 fm2 = internal_force(element, qm2);
        const Vec24 fm1 = internal_force(element, qm1);
        const Vec24 fp1 = internal_force(element, qp1);
        const Vec24 fp2 = internal_force(element, qp2);

        numerical.col(column) =
            (-fp2 + fem::Precision(8) * fp1
                  - fem::Precision(8) * fm1 + fm2)
            / (fem::Precision(12) * step);
    }

    return numerical;
}

fem::Precision relative_matrix_error(const Mat24& a, const Mat24& b) {
    const fem::Precision denominator = std::max({
        a.norm(), b.norm(), std::numeric_limits<fem::Precision>::min()
    });
    return (a - b).norm() / denominator;
}

fem::Precision relative_asymmetry(const Mat24& matrix) {
    const fem::Precision denominator = std::max(
        matrix.norm(), std::numeric_limits<fem::Precision>::min()
    );
    return (matrix - matrix.transpose()).norm() / denominator;
}

TangentCheckResult sweep_tangent(
    fem::model::FRTShellS4& element,
    const Vec24& q,
    const std::string& state_name
) {
    const Mat24 analytic = analytic_tangent(element, q);
    const std::array<fem::Precision, 7> step_sizes{
        1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5
    };

    std::cout << "\nFRT tangent sweep: " << state_name << "\n";
    std::cout << "analytic symmetry error = "
              << std::scientific << std::setprecision(6)
              << relative_asymmetry(analytic) << "\n";
    std::cout << "         h      rel_error      max_abs_err      fd_asym\n";
    std::cout << "--------------------------------------------------------\n";

    TangentCheckResult best;

    for (const fem::Precision h : step_sizes) {
        const Mat24 numerical = finite_difference_tangent_5point(element, q, h);
        const fem::Precision relative_error = relative_matrix_error(analytic, numerical);
        const fem::Precision max_abs_error = (analytic - numerical).cwiseAbs().maxCoeff();
        const fem::Precision fd_asymmetry = relative_asymmetry(numerical);

        std::cout << std::scientific << std::setprecision(6)
                  << std::setw(10) << h
                  << std::setw(15) << relative_error
                  << std::setw(17) << max_abs_error
                  << std::setw(15) << fd_asymmetry << "\n";

        if (relative_error < best.relative_error) {
            best.h = h;
            best.relative_error = relative_error;
            best.max_abs_error = max_abs_error;
            best.fd_asymmetry = fd_asymmetry;
            best.numerical = numerical;
        }
    }

    std::cout << "best h = " << std::scientific << std::setprecision(6) << best.h
              << ", rel_error = " << best.relative_error
              << ", fd_asym = " << best.fd_asymmetry << "\n";
    std::cout << " col node dof       rel_col_err       max_abs_err\n";

    for (fem::Index column = 0; column < NumDofs; ++column) {
        const auto a = analytic.col(column);
        const auto n = best.numerical.col(column);
        const fem::Precision denominator = std::max({
            a.norm(), n.norm(), std::numeric_limits<fem::Precision>::min()
        });
        const fem::Precision relative_column_error = (a - n).norm() / denominator;
        const fem::Precision max_abs_column_error = (a - n).cwiseAbs().maxCoeff();

        std::cout << std::setw(4) << column
                  << std::setw(5) << column / DofsPerNode
                  << std::setw(4) << column % DofsPerNode
                  << std::scientific << std::setprecision(6)
                  << std::setw(18) << relative_column_error
                  << std::setw(18) << max_abs_column_error << "\n";
    }

    return best;
}

Vec24 moderate_state() {
    Vec24 q;
    q <<
         0.00,  0.00,  0.00,   0.05, -0.04,  0.03,
         0.08, -0.03,  0.10,   0.18,  0.07, -0.12,
         0.12,  0.06,  0.18,   0.24, -0.16,  0.11,
        -0.04,  0.05,  0.07,  -0.10,  0.13,  0.20;
    return q;
}

Vec24 large_rotation_state() {
    Vec24 q;
    q <<
         0.00,  0.00,  0.00,   0.40, -0.25,  0.15,
         0.16, -0.08,  0.22,   0.85,  0.35, -0.45,
         0.25,  0.14,  0.35,   1.05, -0.70,  0.55,
        -0.10,  0.12,  0.18,  -0.55,  0.65,  0.90;
    return q;
}

Vec24 mixed_bending_state() {
    Vec24 q;
    q <<
         0.00,  0.00,  0.00,   0.00,  0.00,  0.00,
         0.04,  0.01,  0.20,   0.10, -0.30,  0.08,
         0.06,  0.05,  0.42,   0.35, -0.45,  0.20,
        -0.02,  0.03,  0.18,   0.28, -0.12, -0.18;
    return q;
}

void expect_consistent_tangent(
    fem::model::FRTShellS4& element,
    const Vec24& q,
    const std::string& state_name
) {
    const TangentCheckResult result = sweep_tangent(element, q, state_name);

    EXPECT_LT(result.relative_error, 1e-5)
        << "FRT shell tangent is inconsistent in state " << state_name
        << " (best h=" << result.h << ")";

    EXPECT_LT(result.fd_asymmetry, 1e-5)
        << "Numerical d(f_int)/dq is materially asymmetric in state " << state_name;
}

} // namespace

TEST(FRTShellTangent, Full24x24ModerateState) {
    auto model = build_frt_s4_model();
    auto* element = model._data->elements[0]->as<fem::model::FRTShellS4>();
    ASSERT_NE(element, nullptr);
    expect_consistent_tangent(*element, moderate_state(), "moderate");
}

TEST(FRTShellTangent, Full24x24LargeRotationState) {
    auto model = build_frt_s4_model();
    auto* element = model._data->elements[0]->as<fem::model::FRTShellS4>();
    ASSERT_NE(element, nullptr);
    expect_consistent_tangent(*element, large_rotation_state(), "large_rotation");
}

TEST(FRTShellTangent, Full24x24MixedBendingState) {
    auto model = build_frt_s4_model();
    auto* element = model._data->elements[0]->as<fem::model::FRTShellS4>();
    ASSERT_NE(element, nullptr);
    expect_consistent_tangent(*element, mixed_bending_state(), "mixed_bending");
}
