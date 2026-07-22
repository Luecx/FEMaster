/**
 * @file section_shell.cpp
 * @brief Implements common shell section construction and output-basis handling.
 *
 * The abstract base class contains no constitutive shell law. This translation
 * unit implements only validation, common output transformation, basis
 * construction and diagnostic output shared by every concrete section.
 *
 * @see ShellSection
 *
 * @author Finn Eggers
 * @date 22.07.2026
 */

#include "section_shell.h"

#include "../core/logging.h"

#include <Eigen/Geometry>

#include <sstream>
#include <utility>

namespace fem {

ShellSection::ShellSection(
    material::Material::Ptr    material,
    model::ElementRegion::Ptr  region,
    Precision                  thickness,
    cos::CoordinateSystem::Ptr orientation,
    Index                      csys_axis
)
    : thickness_  (thickness),
      orientation_(std::move(orientation)),
      csys_axis_  (csys_axis) {
    // A non-positive thickness invalidates constitutive integration, mass
    // integration and every physical thickness-coordinate mapping.
    logging::error(thickness_ > Precision(0),
        "ShellSection: thickness must be positive");

    // The external one-based CSYSAXIS value has already been converted to this
    // zero-based index before construction.
    logging::error(csys_axis_ >= 0 && csys_axis_ < 3,
        "ShellSection: CSYSAXIS must be 1, 2 or 3");

    // Store the associations inherited from the generic Section base class.
    material_ = std::move(material);
    region_   = std::move(region);
}

ShellStressResultants ShellSection::evaluate_output_resultants(
    const Vec3&                   position_reference,
    const Mat3&                   shell_basis_global,
    const ShellGeneralizedStrain& strain_shell,
    bool                          use_green_lagrange
) const {
    // Every concrete section returns the element-facing response in the
    // geometric shell basis. Keeping this call in the base class avoids a
    // duplicate resultant-output implementation in every section type.
    ShellStressResultants resultants_shell;
    Mat8                  tangent_shell;

    evaluate(
        position_reference,
        shell_basis_global,
        strain_shell,
        use_green_lagrange,
        resultants_shell,
        tangent_shell
    );

    // Construct the basis used for the stored membrane-force, moment and
    // transverse-shear components.
    const Mat3 result_basis_global =
        stress_resultant_basis(position_reference, shell_basis_global);

    // Express the two result axes in geometric shell coordinates:
    //
    //     R = Q_shell^T Q_result.
    //
    // The columns of R are the target axes written in the current basis.
    const Mat2 result_axes_in_shell =
        shell_basis_global.template block<3, 2>(0, 0).transpose()
        * result_basis_global.template block<3, 2>(0, 0);

    // Rotate the physical resultant components into the configured output basis.
    return resultants_shell.transformed(result_axes_in_shell);
}

Mat3 ShellSection::stress_basis(
    const Vec3& position_reference,
    const Mat3& shell_basis_global
) const {
    // Without a user orientation, physical stresses are ordinary global tensor
    // components and therefore require no local basis construction.
    if (!orientation_) {
        return Mat3::Identity();
    }

    const Precision projection_tolerance = Precision(1e-6);

    // Normalize defensively so the projection tolerance has a clear geometric
    // meaning independent of small numerical basis errors.
    const Vec3 normal = shell_basis_global.col(2).normalized();

    // Curvilinear coordinate systems depend on position. Convert the physical
    // reference point first and then evaluate the local coordinate-system axes.
    const Vec3 point_local      = orientation_->to_local(position_reference);
    const Mat3 orientation_axes = orientation_->get_axes(point_local);
    const Vec3 source_axis      = orientation_axes.col(csys_axis_);

    // Remove the normal component of the selected axis:
    //
    //     e1_raw = a - (a . n) n.
    Vec3 stress_e1 = source_axis - normal * source_axis.dot(normal);

    // An explicitly selected axis must be valid everywhere on the section. A
    // silent fallback would make the requested material/output orientation
    // ineffective and is therefore forbidden.
    logging::error(stress_e1.norm() > projection_tolerance,
        "ShellSection: CSYSAXIS ", csys_axis_ + 1,
        " is perpendicular to the shell surface and cannot define the stress basis");

    stress_e1.normalize();

    // Complete the right-handed orthonormal basis with the shell normal.
    const Vec3 stress_e2 = normal.cross(stress_e1).normalized();

    Mat3 basis;
    basis.col(0) = stress_e1;
    basis.col(1) = stress_e2;
    basis.col(2) = normal;

    return basis;
}

Mat3 ShellSection::stress_resultant_basis(
    const Vec3& position_reference,
    const Mat3& shell_basis_global
) const {
    // One explicit orientation defines both physical stress and generalized
    // stress-resultant components.
    if (orientation_) {
        return stress_basis(position_reference, shell_basis_global);
    }

    const Precision projection_tolerance = Precision(1e-6);
    const Vec3      normal               = shell_basis_global.col(2).normalized();

    // Resultants are stored in local membrane, moment and transverse-shear
    // blocks. A deterministic tangent basis is therefore required before
    // neighboring element contributions can be averaged component-wise.
    Vec3 source_axis  = Vec3::UnitX();
    Vec3 resultant_e1 = source_axis - normal * source_axis.dot(normal);

    // Global Y is only an automatic output fallback. Explicit user orientations
    // never fall back silently.
    if (resultant_e1.norm() <= projection_tolerance) {
        source_axis  = Vec3::UnitY();
        resultant_e1 = source_axis - normal * source_axis.dot(normal);
    }

    // Valid three-dimensional shell geometry cannot make both global X and Y
    // parallel to the same normal. This check catches invalid/non-finite bases.
    logging::error(resultant_e1.norm() > projection_tolerance,
        "ShellSection: global X and Y cannot define a tangential stress-resultant basis");

    resultant_e1.normalize();

    const Vec3 resultant_e2 = normal.cross(resultant_e1).normalized();

    Mat3 basis;
    basis.col(0) = resultant_e1;
    basis.col(1) = resultant_e2;
    basis.col(2) = normal;

    return basis;
}

void ShellSection::info() {
    // All concrete sections share these model associations and orientation data.
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material   : ", (material_    ? material_    ->name : "-"));
    logging::info(true, "   Region     : ", (region_      ? region_      ->name : "-"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
    logging::info(orientation_ != nullptr, "   CSYSAXIS   : ", csys_axis_ + 1);
    logging::info(true, "   Thickness  : ", thickness_);
}

std::string ShellSection::str() const {
    std::ostringstream os;

    // Keep the compact representation stable for logs and model summaries.
    os << "ShellSection: t=" << thickness_
       << ", material="      << (material_    ? material_    ->name : std::string("-"))
       << ", orientation="   << (orientation_ ? orientation_->name : std::string("-"))
       << ", csysaxis="      << (orientation_ ? std::to_string(csys_axis_ + 1) : std::string("-"))
       << ", region="        << (region_      ? region_      ->name : std::string("-"))
       << " ("               << (region_      ? static_cast<int>(region_->size()) : 0) << ")";

    return os.str();
}

} // namespace fem
