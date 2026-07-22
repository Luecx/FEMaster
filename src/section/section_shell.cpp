/**
 * @file section_shell.cpp
 * @brief Implements common shell section behavior and basis transformations.
 *
 * @see src/section/section_shell.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "section_shell.h"

#include "../core/logging.h"

#include <Eigen/Geometry>

#include <sstream>
#include <utility>

namespace fem {

ShellSection::ShellSection(material::Material::Ptr    material,
                           model::ElementRegion::Ptr  region,
                           Precision                  thickness,
                           cos::CoordinateSystem::Ptr orientation,
                           Index                      csys_axis)
    : thickness_  (thickness),
      orientation_(std::move(orientation)),
      csys_axis_  (csys_axis) {
    logging::error(thickness_ > Precision(0), "ShellSection: thickness must be positive");
    logging::error(csys_axis_ >= 0 && csys_axis_ < 3,
        "ShellSection: CSYSAXIS must be 1, 2 or 3");

    material_ = std::move(material);
    region_   = std::move(region);
}

void ShellSection::evaluate(const Vec3&                   position_reference,
                            const Mat3&                   shell_basis_global,
                            const ShellGeneralizedStrain& strain,
                            bool                          use_green_lagrange,
                            ShellStressResultants&        resultants,
                            Mat8&                         tangent) const {
    if (!orientation_) {
        evaluate_material(strain, use_green_lagrange, resultants, tangent);
        return;
    }

    // Build the material basis in the supplied shell plane
    const Mat3 material_basis_global = output_basis_global(position_reference, shell_basis_global);

    // Transform generalized shell strains into the material basis
    Mat2 material_to_shell;
    material_to_shell(0, 0) = shell_basis_global.col(0).dot(material_basis_global.col(0));
    material_to_shell(1, 0) = shell_basis_global.col(1).dot(material_basis_global.col(0));
    material_to_shell(0, 1) = shell_basis_global.col(0).dot(material_basis_global.col(1));
    material_to_shell(1, 1) = shell_basis_global.col(1).dot(material_basis_global.col(1));

    const Mat8                   strain_transform = ShellGeneralizedStrain::transformation(material_to_shell);
    const ShellGeneralizedStrain strain_material(strain_transform * strain.values());

    ShellStressResultants resultants_material;
    Mat8                  tangent_material;
    evaluate_material(
        strain_material,
        use_green_lagrange,
        resultants_material,
        tangent_material
    );

    // Return resultants and tangent in the supplied shell basis. Stress
    // resultants transform in the inverse direction with the physical-shear
    // convention, which is dual to the generalized strain transformation.
    const Mat2 shell_to_material    = material_to_shell.transpose();
    const Mat8 resultants_transform = ShellStressResultants::transformation(shell_to_material);

    resultants = ShellStressResultants(resultants_transform * resultants_material.values());
    tangent    = resultants_transform * tangent_material * strain_transform;
}

/**
 * Evaluates generalized shell resultants in the output basis.
 *
 * With an assigned orientation, the generalized strain is transformed into the
 * projected material basis and the section is evaluated directly there. Without
 * an orientation, the section is evaluated in the local shell basis and only the
 * output resultants are transformed into a stable in-plane basis constructed by
 * projecting global X into the shell plane, falling back to global Y when X is
 * nearly normal to the shell.
 *
 * @param position_reference Physical reference position of the evaluation point.
 * @param shell_basis_global Orthonormal shell basis defining the input strain.
 * @param strain Generalized shell strain in the supplied shell basis.
 * @param use_green_lagrange Select finite-strain material evaluation.
 * @return Generalized shell resultants in the output basis.
 */
ShellStressResultants ShellSection::compute_resultants(const Vec3&                   position_reference,
                                                       const Mat3&                   shell_basis_global,
                                                       const ShellGeneralizedStrain& strain,
                                                       bool                          use_green_lagrange) const {
    ShellStressResultants resultants_output;
    Mat8                  tangent_output;

    if (orientation_) {
        // Evaluate directly in the material/output basis when the user supplied
        // an orientation.
        const Mat2 output_to_shell  = output_to_shell_rotation(position_reference, shell_basis_global);
        const Mat8 strain_transform = ShellGeneralizedStrain::transformation(output_to_shell);
        const ShellGeneralizedStrain strain_output(strain_transform * strain.values());

        evaluate_material(
            strain_output,
            use_green_lagrange,
            resultants_output,
            tangent_output
        );

        return resultants_output;
    }

    ShellStressResultants resultants_shell;

    // Without an orientation the material model is evaluated in the element's
    // local shell basis. Only the output resultants are rotated into the stable
    // global-axis projection basis.
    evaluate_material(
        strain,
        use_green_lagrange,
        resultants_shell,
        tangent_output
    );

    const Mat2 output_to_shell = output_to_shell_rotation(position_reference, shell_basis_global);
    return resultants_shell.transformed(output_to_shell);
}

/**
 * Reports that direct material stress recovery is unavailable for this section.
 *
 * Sections without through-thickness material integration, such as prescribed
 * ABD sections, do not contain enough constitutive information to evaluate a
 * material-point stress through the material model.
 *
 * @return Never returns during normal execution because the section cannot
 *         compute material stresses.
 */
VolumeStressCauchy ShellSection::compute_stress(const Vec3&                   position_reference,
                                                const Mat3&                   shell_basis_global,
                                                const ShellGeneralizedStrain& strain,
                                                Precision                     z,
                                                bool                          use_green_lagrange,
                                                Mat3                          deformation_gradient) const {
    (void) position_reference;
    (void) shell_basis_global;
    (void) strain;
    (void) z;
    (void) use_green_lagrange;
    (void) deformation_gradient;

    logging::error(false,
        "ShellSection: material stress recovery requires a through-thickness integrated section");

    return VolumeStressCauchy();
}

/**
 * Builds the global basis used for material orientation and resultant output.
 *
 * If a coordinate system is assigned, the selected `CSYSAXIS` is projected into
 * the shell plane. A nearly normal selected axis is treated as an invalid model
 * definition because it cannot define a stable material direction.
 *
 * Without a coordinate system, the output basis is made deterministic by first
 * projecting global X into the shell plane and then falling back to global Y
 * when X is nearly perpendicular to the shell.
 *
 * @param position_reference Physical reference position of the evaluation point.
 * @param shell_basis_global Orthonormal shell basis whose third column is the shell normal.
 * @return Orthonormal global output basis.
 */
Mat3 ShellSection::output_basis_global(const Vec3& position_reference,
                                       const Mat3& shell_basis_global) const {
    const Precision projection_tolerance = Precision(1e-6);
    const Vec3      normal               = shell_basis_global.col(2).normalized();

    // Select the source direction that defines the first output axis
    Vec3 source_axis;
    if (orientation_) {
        const Vec3 point_local      = orientation_->to_local(position_reference);
        const Mat3 orientation_axes = orientation_->get_axes(point_local);
        source_axis                 = orientation_axes.col(csys_axis_);
    } else {
        source_axis << Precision(1), Precision(0), Precision(0);
    }

    // Project the selected direction into the tangent plane
    Vec3 output_e1 = source_axis - normal * source_axis.dot(normal);

    if (!orientation_ && output_e1.norm() <= projection_tolerance) {
        source_axis << Precision(0), Precision(1), Precision(0);
        output_e1 = source_axis - normal * source_axis.dot(normal);
    }

    if (orientation_) {
        logging::error(output_e1.norm() > projection_tolerance,
            "ShellSection: CSYSAXIS ", csys_axis_ + 1,
            " cannot be projected onto some shell elements because it is perpendicular to the shell surface");
    } else {
        logging::error(output_e1.norm() > projection_tolerance,
            "ShellSection: global X and Y axes cannot be projected onto some shell elements");
    }

    output_e1.normalize();

    const Vec3 output_e2 = normal.cross(output_e1).normalized();

    Mat3 output_basis;
    output_basis.col(0) = output_e1;
    output_basis.col(1) = output_e2;
    output_basis.col(2) = normal;

    return output_basis;
}

/**
 * Builds the two-dimensional rotation from the local shell basis to the output
 * basis.
 *
 * Columns contain the output in-plane basis vectors expressed in the supplied
 * shell basis, matching the convention used by the shell resultant
 * transformation matrices.
 *
 * @param position_reference Physical reference position of the evaluation point.
 * @param shell_basis_global Orthonormal shell basis used by the shell element.
 * @return Output basis vectors expressed in the shell basis.
 */
Mat2 ShellSection::output_to_shell_rotation(const Vec3& position_reference,
                                            const Mat3& shell_basis_global) const {
    const Mat3 output_basis = output_basis_global(position_reference, shell_basis_global);

    Mat2 output_to_shell;
    output_to_shell(0, 0) = shell_basis_global.col(0).dot(output_basis.col(0));
    output_to_shell(1, 0) = shell_basis_global.col(1).dot(output_basis.col(0));
    output_to_shell(0, 1) = shell_basis_global.col(0).dot(output_basis.col(1));
    output_to_shell(1, 1) = shell_basis_global.col(1).dot(output_basis.col(1));

    return output_to_shell;
}

void ShellSection::info() {
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material   : ", (material_    ? material_   ->name : "-"));
    logging::info(true, "   Region     : ", (region_      ? region_     ->name : "-"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
    logging::info(orientation_ != nullptr, "   CSYSAXIS   : ", csys_axis_ + 1);
    logging::info(true, "   Thickness  : ", thickness_);
}

std::string ShellSection::str() const {
    std::ostringstream os;

    os << "ShellSection: t=" << thickness_
       << ", material="      << (material_     ? material_   ->name : std::string("-"))
       << ", orientation="   << (orientation_  ? orientation_->name : std::string("-"))
       << ", csysaxis="      << (orientation_  ? std::to_string(csys_axis_ + 1) : std::string("-"))
       << ", region="        << (region_       ? region_     ->name : std::string("-"))
       << " ("               << (region_       ? static_cast<int>(region_->size()) : 0) << ")";

    return os.str();
}
} // namespace fem
