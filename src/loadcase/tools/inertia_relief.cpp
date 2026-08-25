/**
 * @file inertia_relief.cpp
 * @brief Implements inertia relief for unconstrained structural systems.
 *
 * Inertia relief replaces missing rigid-body constraints by equivalent
 * translational and rotational inertia loads. The rigid-body accelerations are
 * selected such that the externally applied resultant force and moment are
 * balanced by the corresponding inertia forces.
 *
 * The implementation performs the following operations:
 *
 * 1. Integrate the compiled structural mass and first mass moment.
 * 2. Optionally include auxiliary PointElements created by native POINTMASS.
 * 3. Compute the center of gravity.
 * 4. Reduce the external nodal loads to a resultant force and moment.
 * 5. Integrate the inertia tensor about the center of gravity.
 * 6. Optionally include translational and rotary inertia of auxiliary
 *    PointElements.
 * 7. Solve the translational and rotational rigid-body equations.
 * 8. Assemble the resulting inertia loads through `bc::InertialLoad`.
 * 9. Verify the final rigid-body equilibrium.
 *
 * Regular compiled PointElements are structural elements and therefore
 * participate in the normal element integration path. Auxiliary PointElements
 * are stored separately in `ModelData::point_elements`; they participate only
 * when `consider_point_masses` is enabled, preserving the existing FEMaster
 * inertia-relief option for native POINTMASS definitions.
 *
 * The rotational equation may be singular when the mass distribution does not
 * provide inertia about one or more axes. The symmetric inertia tensor is
 * therefore diagonalized, and a pseudo-inverse is formed from its principal
 * inertia values.
 *
 * @see apply_inertia_relief
 * @see bc::InertialLoad
 * @see model::PointElement
 * @see PointMassSection
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "inertia_relief.h"

#include "../../bc/load_inertial.h"
#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/element/point.h"
#include "../../section/section_point_mass.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <utility>

namespace fem {
namespace {

/**
 * Constructs the unit-mass inertia tensor associated with a point.
 *
 * Let `r` denote the position of a unit point mass relative to the selected
 * reference point. Its inertia tensor is
 *
 *     J(r) = (r^T r) I - r r^T.
 *
 * In Cartesian coordinates,
 *
 *           [ y² + z²    -xy        -xz     ]
 *     J  =  [   -xy     x² + z²     -yz     ].
 *           [   -xz       -yz      x² + y²  ]
 *
 * Multiplying this tensor by a point mass gives its physical translational
 * inertia contribution. Integrating it with respect to the structural mass
 * measure gives the inertia tensor of a distributed body.
 *
 * @param relative_position Position relative to the inertia reference point.
 *
 * @return Symmetric inertia tensor for a unit point mass.
 */
Mat3 inertia_about_point(const Vec3& relative_position) {
    // Extract the Cartesian lever-arm components
    const Precision x = relative_position(0);
    const Precision y = relative_position(1);
    const Precision z = relative_position(2);

    // Precompute all quadratic products required by
    //
    //     J = (r^T r) I - r r^T
    //
    // so the tensor can be assembled without temporary matrix products
    const Precision xx = x * x;
    const Precision yy = y * y;
    const Precision zz = z * z;
    const Precision xy = x * y;
    const Precision xz = x * z;
    const Precision yz = y * z;

    // Assemble the symmetric unit-mass inertia tensor explicitly
    Mat3 inertia = Mat3::Zero();

    inertia(0, 0) =  yy + zz;
    inertia(1, 1) =  xx + zz;
    inertia(2, 2) =  xx + yy;

    inertia(0, 1) = -xy;
    inertia(1, 0) = -xy;

    inertia(0, 2) = -xz;
    inertia(2, 0) = -xz;

    inertia(1, 2) = -yz;
    inertia(2, 1) = -yz;

    return inertia;
}

/**
 * Integrates the compiled structural mass and first mass moment.
 *
 * The structural integration routines use the same physical mass measure as
 * the inertial-load implementation. Therefore,
 *
 *     m     = integral dm,
 *     s_m   = integral x dm.
 *
 * This includes regular compiled PointElements because their density-scaled
 * integration interprets `PointMassSection::mass_` as a concentrated mass
 * measure. Auxiliary PointElements are not part of `ModelData::elements` and
 * are added separately when requested.
 *
 * @param model_data Global model data containing the compiled element collection.
 * @param mass Accumulated structural mass.
 * @param first_mass_moment Accumulated first mass moment about the global origin.
 */
void accumulate_structural_mass(const model::ModelData& model_data,
                                Precision&              mass,
                                Vec3&                   first_mass_moment) {
    // Integrate all regular compiled structural elements
    for (const auto& element : model_data.elements) {
        if (!element) continue;

        auto* structural_element = element->as<model::StructuralElement>();
        if (!structural_element) continue;

        // Integrating one gives the physical mass contribution of the element.
        mass += structural_element->integrate_scalar_field(true,
            [](const Vec3&) { return Precision(1); });

        // Integrating position gives the first mass moment about the origin.
        first_mass_moment += structural_element->integrate_vector_field(true,
            [](const Vec3& position) { return position; });
    }
}

/**
 * Adds the auxiliary native point masses to total mass and first mass moment.
 *
 * Native `POINTMASS, NSET=...` definitions are represented by one auxiliary
 * `PointElement` per target node and stored in `ModelData::point_elements` after
 * compilation. They are deliberately excluded from the regular element array,
 * so `consider_point_masses` adds them explicitly here.
 *
 * For every auxiliary point element with mass `m_i` at `x_i`,
 *
 *     m   += m_i,
 *     s_m += m_i x_i.
 *
 * Rotary inertia does not affect the center of gravity.
 *
 * @param model_data Global model data containing auxiliary point elements.
 * @param positions Global nodal coordinate field.
 * @param mass Accumulated total mass.
 * @param first_mass_moment Accumulated first mass moment.
 */
void accumulate_point_mass_mass(const model::ModelData& model_data,
                                const model::Field&     positions,
                                Precision&              mass,
                                Vec3&                   first_mass_moment) {
    // `consider_point_masses` refers exactly to these post-compile auxiliary
    // PointElements. Do not inspect the generic feature collection here.
    for (const auto& element : model_data.point_elements) {
        if (!element) continue;

        const auto* point = element->as<model::PointElement>();
        logging::error(point != nullptr,
            "InertiaRelief: auxiliary point-element storage contains a non-point element");
        logging::error(point->_section != nullptr,
            "InertiaRelief: auxiliary point element ", point->elem_id, " has no section");

        const auto* section = point->_section->as<PointMassSection>();
        logging::error(section != nullptr,
            "InertiaRelief: auxiliary point element ", point->elem_id,
            " has a non-point-mass section");

        if (section->mass_ == Precision(0)) continue;

        const ID node_id = point->nodes()[0];
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
            "InertiaRelief: point-mass node ", node_id,
            " is outside the positions field with ", positions.rows, " rows");

        const Index node     = static_cast<Index>(node_id);
        const Vec3  position = positions.row_vec3(node);

        mass              += section->mass_;
        first_mass_moment += section->mass_ * position;
    }
}

/**
 * Integrates the inertia tensor of regular compiled structural elements about
 * the center of gravity.
 *
 * For every mass point with relative position
 *
 *     r = x - x_c,
 *
 * the translational inertia contribution is
 *
 *     dI = [(r^T r) I - r r^T] dm.
 *
 * Regular PointElements contribute through the same mass-weighted integration.
 * Their diagonal rotary inertia is added explicitly because it is an intrinsic
 * nodal inertia and is not part of the geometric parallel-axis kernel.
 *
 * @param model_data Global model data containing compiled structural elements.
 * @param center_of_gravity Global center of gravity.
 * @param inertia Accumulated inertia tensor about the center of gravity.
 */
void accumulate_structural_inertia(const model::ModelData& model_data,
                                   const Vec3&             center_of_gravity,
                                   Mat3&                   inertia) {
    // Integrate the translational inertia of every compiled structural element
    for (const auto& element : model_data.elements) {
        if (!element) continue;

        auto* structural_element = element->as<model::StructuralElement>();
        if (!structural_element) continue;

        inertia += structural_element->integrate_tensor_field(true,
            [&](const Vec3& position) {
                const Vec3 relative_position = position - center_of_gravity;
                return inertia_about_point(relative_position);
            });

        // PointMassSection may additionally carry intrinsic rotary inertia. It
        // contributes directly about the nodal axes and must not be multiplied
        // by the point position or translational mass again.
        if (const auto* point = element->as<model::PointElement>()) {
            if (!point->_section) continue;

            const auto* section = point->_section->as<PointMassSection>();
            logging::error(section != nullptr,
                "InertiaRelief: point element ", point->elem_id,
                " has a non-point-mass section");
            inertia.diagonal() += section->rotary_inertia_;
        }
    }
}

/**
 * Adds the inertia of auxiliary native PointElements about the center of gravity.
 *
 * Each auxiliary element contributes
 *
 *     I_trans = m [(r^T r) I - r r^T]
 *
 * from its concentrated translational mass and
 *
 *     I_rot = diag(Jx, Jy, Jz)
 *
 * from `PointMassSection::rotary_inertia_`. These objects are considered only
 * when the caller enables `consider_point_masses`.
 *
 * @param model_data Global model data containing auxiliary point elements.
 * @param positions Global nodal coordinate field.
 * @param center_of_gravity Global center of gravity.
 * @param inertia Accumulated inertia tensor.
 */
void accumulate_point_mass_inertia(const model::ModelData& model_data,
                                   const model::Field&     positions,
                                   const Vec3&             center_of_gravity,
                                   Mat3&                   inertia) {
    // Traverse the same auxiliary PointElements used by InertialLoad when
    // `consider_point_masses_` is enabled.
    for (const auto& element : model_data.point_elements) {
        if (!element) continue;

        const auto* point = element->as<model::PointElement>();
        logging::error(point != nullptr,
            "InertiaRelief: auxiliary point-element storage contains a non-point element");
        logging::error(point->_section != nullptr,
            "InertiaRelief: auxiliary point element ", point->elem_id, " has no section");

        const auto* section = point->_section->as<PointMassSection>();
        logging::error(section != nullptr,
            "InertiaRelief: auxiliary point element ", point->elem_id,
            " has a non-point-mass section");

        const ID node_id = point->nodes()[0];
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
            "InertiaRelief: point-mass node ", node_id,
            " is outside the positions field with ", positions.rows, " rows");

        const Index node              = static_cast<Index>(node_id);
        const Vec3  position          = positions.row_vec3(node);
        const Vec3  relative_position = position - center_of_gravity;

        inertia += section->mass_ * inertia_about_point(relative_position);
        inertia.diagonal() += section->rotary_inertia_;
    }
}

} // namespace

/**
 * Applies inertia relief to an assembled nodal load field.
 *
 * The external nodal loads are reduced to a force resultant `F` and a moment
 * resultant `M_c` about the center of gravity. The rigid-body accelerations are
 * then determined from
 *
 *     m a       = F,
 *     I_c alpha = M_c.
 *
 * The translational equation is solved directly. The rotational equation is
 * solved through the pseudo-inverse of the symmetric inertia tensor because
 * geometrically degenerate models may have one or more vanishing principal
 * inertias.
 *
 * Let
 *
 *     I_c = Q diag(I_1, I_2, I_3) Q^T
 *
 * denote the principal-axis decomposition. The pseudo-inverse is constructed as
 *
 *     I_c^+ = Q diag(I_1^+, I_2^+, I_3^+) Q^T,
 *
 * with
 *
 *                 1 / I_i,  |I_i| > epsilon_I,
 *     I_i^+ = {
 *                 0,        otherwise.
 *
 * Consequently, no angular acceleration is generated about a direction for
 * which the model possesses no meaningful rotational inertia.
 *
 * Compiled structural elements always participate. Auxiliary PointElements from
 * native POINTMASS commands participate in the mass properties and in the
 * resulting `bc::InertialLoad` only when `consider_point_masses` is enabled.
 *
 * @param model_data Global model data.
 * @param global_load_mat Nodal load field modified directly by inertia relief.
 * @param consider_point_masses Include auxiliary native PointElements and their
 *                              translational/rotary inertia.
 */
void apply_inertia_relief(model::ModelData& model_data,
                          model::Field&     global_load_mat,
                          bool              consider_point_masses) {
    // Absolute lower bound used when constructing the pseudo-inverse of the
    // rotational inertia tensor
    constexpr Precision minimum_inertia_eigenvalue = Precision(1e-12);

    // Relative and absolute tolerances used for the final rigid-body
    // equilibrium verification
    constexpr Precision relative_balance_tolerance = Precision(1e-8);
    constexpr Precision absolute_balance_tolerance = Precision(1e-10);

    // Formatting parameters for the diagnostic table
    constexpr int log_width     = 12;
    constexpr int log_precision = 3;

    // ---------------------------------------------------------------------
    // Validate the model fields required by all subsequent calculations
    // ---------------------------------------------------------------------

    logging::error(model_data.positions != nullptr,
        "InertiaRelief: positions not initialized");

    auto& positions = *model_data.positions;

    logging::error(global_load_mat.rows == positions.rows,
        "InertiaRelief: global_load_mat.rows != positions.rows (",
        global_load_mat.rows, " vs ", positions.rows, ")");
    logging::error(global_load_mat.components >= 3,
        "InertiaRelief: global_load_mat.components must be at least 3");

    /**
     * Reduces a nodal load field to its rigid-body resultants.
     *
     * The total force is
     *
     *     F = sum_i f_i.
     *
     * The total moment about a reference point `c` is
     *
     *     M_c = sum_i [m_i + (x_i - c) x f_i],
     *
     * where `m_i` denotes an optional directly applied nodal moment. Nodal
     * moments are included only when the field provides at least six
     * components.
     */
    const auto compute_resultants =
        [&](const model::Field& loads, const Vec3& reference_point)
        -> std::pair<Vec3, Vec3> {
            Vec3 resultant_force  = Vec3::Zero();
            Vec3 resultant_moment = Vec3::Zero();

            const bool has_nodal_moments = loads.components >= 6;

            // Sum the contribution of every node
            for (Index node = 0; node < loads.rows; ++node) {
                // Components zero to two represent the nodal force vector
                const Vec3 nodal_force{
                    loads(node, 0),
                    loads(node, 1),
                    loads(node, 2)
                };

                // Components three to five optionally represent a directly
                // applied nodal moment
                Vec3 nodal_moment = Vec3::Zero();

                if (has_nodal_moments) {
                    nodal_moment = Vec3{
                        loads(node, 3),
                        loads(node, 4),
                        loads(node, 5)
                    };
                }

                // The lever arm must be measured from the same reference point
                // about which the rotational equation will be formulated
                const Vec3 position          = positions.row_vec3(node);
                const Vec3 relative_position = position - reference_point;

                // Accumulate
                //
                //     F   <- F   + f_i,
                //     M_c <- M_c + m_i + r_i x f_i
                resultant_force  += nodal_force;
                resultant_moment += nodal_moment + relative_position.cross(nodal_force);
            }

            return {resultant_force, resultant_moment};
        };

    // ---------------------------------------------------------------------
    // Integrate mass and determine the global center of gravity
    // ---------------------------------------------------------------------

    Precision total_mass        = Precision(0);
    Vec3      first_mass_moment = Vec3::Zero();

    // Regular compiled structural elements, including compiled PointElements,
    // always participate in inertia relief.
    accumulate_structural_mass(model_data, total_mass, first_mass_moment);

    // Native post-compile POINTMASS elements participate only when requested.
    if (consider_point_masses) {
        accumulate_point_mass_mass(model_data, positions, total_mass, first_mass_moment);
    }

    // Validate all integrated mass quantities together before dividing by the
    // total mass
    logging::error(total_mass > Precision(0) && std::isfinite(total_mass),
        "InertiaRelief: total mass is zero or invalid");
    logging::error(first_mass_moment.allFinite(),
        "InertiaRelief: first mass moment contains NaN or Inf");

    // The center of gravity follows from
    //
    //     x_c = (integral x dm) / (integral dm)
    const Vec3 center_of_gravity = first_mass_moment / total_mass;

    logging::error(center_of_gravity.allFinite(),
        "InertiaRelief: center of gravity contains NaN or Inf");

    // ---------------------------------------------------------------------
    // Reduce the external loading to force and moment resultants
    // ---------------------------------------------------------------------

    // The moment is evaluated about the center of gravity. This decouples the
    // rigid-body translational and rotational equations:
    //
    //     m a       = F,
    //     I_c alpha = M_c.
    const auto [external_force, external_moment] =
        compute_resultants(global_load_mat, center_of_gravity);

    logging::error(external_force.allFinite() && external_moment.allFinite(),
        "InertiaRelief: external resultants contain NaN or Inf");

    // ---------------------------------------------------------------------
    // Integrate the inertia tensor about the center of gravity
    // ---------------------------------------------------------------------

    Mat3 inertia = Mat3::Zero();

    // Regular compiled structural contribution, including PointElements.
    accumulate_structural_inertia(model_data, center_of_gravity, inertia);

    // Optional native post-compile PointElement contribution.
    if (consider_point_masses) {
        accumulate_point_mass_inertia(model_data, positions, center_of_gravity, inertia);
    }

    logging::error(inertia.allFinite(),
        "InertiaRelief: inertia tensor contains NaN or Inf");

    // ---------------------------------------------------------------------
    // Solve the translational rigid-body equation
    // ---------------------------------------------------------------------

    // Because the moment is formulated about the center of gravity, rigid-body
    // translation is governed independently by
    //
    //     m a = F
    const Vec3 translational_acceleration = external_force / total_mass;

    // ---------------------------------------------------------------------
    // Diagonalize the rotational inertia tensor
    // ---------------------------------------------------------------------

    // The inertia tensor is real and symmetric by construction. A
    // SelfAdjointEigenSolver therefore provides orthonormal principal axes and
    // real principal inertia values:
    //
    //     I_c = Q Lambda Q^T
    Eigen::SelfAdjointEigenSolver<Mat3> eigen_solver(inertia);

    logging::error(eigen_solver.info() == Eigen::Success,
        "InertiaRelief: eigen decomposition of inertia tensor failed");

    const Vec3 principal_inertias = eigen_solver.eigenvalues();
    const Mat3 principal_axes     = eigen_solver.eigenvectors();

    logging::error(principal_inertias.allFinite() && principal_axes.allFinite(),
        "InertiaRelief: principal inertia decomposition contains NaN or Inf");

    // Use the largest absolute principal inertia as the scale for identifying
    // numerically vanishing rotational modes
    const Precision maximum_principal_inertia =
        principal_inertias.cwiseAbs().maxCoeff();

    if (!std::isfinite(maximum_principal_inertia) ||
        maximum_principal_inertia <= Precision(0)) {
        logging::info(true, inertia);
    }

    logging::error(std::isfinite(maximum_principal_inertia) && maximum_principal_inertia > Precision(0),
        "InertiaRelief: inertia tensor is zero or invalid");

    // Combine an absolute cutoff with a tensor-relative cutoff:
    //
    //     epsilon_I = max(epsilon_abs, 1e-12 max_i |I_i|)
    //
    // This prevents inversion of both exactly zero and numerically negligible
    // principal inertia values.
    const Precision inertia_cutoff = std::max(
        minimum_inertia_eigenvalue,
        Precision(1e-12) * maximum_principal_inertia
    );

    // ---------------------------------------------------------------------
    // Construct the pseudo-inverse in principal coordinates
    // ---------------------------------------------------------------------

    Vec3 inverse_principal_inertias = Vec3::Zero();

    for (Index direction = 0; direction < 3; ++direction) {
        const Precision principal_inertia = principal_inertias(direction);

        // Supported rotational directions receive the normal reciprocal
        //
        //     I_i^+ = 1 / I_i
        if (std::abs(principal_inertia) > inertia_cutoff) {
            inverse_principal_inertias(direction) = Precision(1) / principal_inertia;
        }

        // Unsupported or numerically singular directions retain
        //
        //     I_i^+ = 0
        //
        // and therefore generate no angular acceleration component.
    }

    // Transform the external moment into principal coordinates, divide every
    // supported component by its principal inertia and transform the resulting
    // angular acceleration back into global coordinates:
    //
    //     alpha = Q Lambda^+ Q^T M_c
    const Vec3 angular_acceleration =
        principal_axes *
        inverse_principal_inertias.asDiagonal() *
        principal_axes.transpose() *
        external_moment;

    // Group the final acceleration validation after both rigid-body equations
    // have been solved
    logging::error(translational_acceleration.allFinite() && angular_acceleration.allFinite(),
        "InertiaRelief: computed accelerations contain NaN or Inf");

    // ---------------------------------------------------------------------
    // Report the original resultants and computed accelerations
    // ---------------------------------------------------------------------

    logging::info(true, "InertiaRelief:");
    logging::info(true, "              ",
        std::setw(log_width), "x",
        std::setw(log_width), "y",
        std::setw(log_width), "z",
        std::setw(log_width), "rx",
        std::setw(log_width), "ry",
        std::setw(log_width), "rz");
    logging::info(true, "initial F/M  :", std::setprecision(log_precision),
        std::setw(log_width), external_force(0),
        std::setw(log_width), external_force(1),
        std::setw(log_width), external_force(2),
        std::setw(log_width), external_moment(0),
        std::setw(log_width), external_moment(1),
        std::setw(log_width), external_moment(2));
    logging::info(true, "using a/alpha:", std::setprecision(log_precision),
        std::setw(log_width), translational_acceleration(0),
        std::setw(log_width), translational_acceleration(1),
        std::setw(log_width), translational_acceleration(2),
        std::setw(log_width), angular_acceleration(0),
        std::setw(log_width), angular_acceleration(1),
        std::setw(log_width), angular_acceleration(2));

    // ---------------------------------------------------------------------
    // Assemble and apply the equivalent inertia loads
    // ---------------------------------------------------------------------

    // InertialLoad assembles both translational and rotational inertia
    // contributions into a six-component nodal load field
    model::Field inertial_loads =
        model_data.create_field_("_TEMP", model::FieldDomain::NODE, 6);

    inertial_loads.set_zero();

    bc::InertialLoad inertia_load;

    // The accelerations are defined at the center of gravity, which is also the
    // reference point used for the external moment resultant
    inertia_load.center_                = center_of_gravity;
    inertia_load.center_acc_            = translational_acceleration;
    inertia_load.alpha_                 = angular_acceleration;

    // Inertia relief currently considers no prescribed rigid-body angular
    // velocity and therefore produces no centrifugal contribution
    inertia_load.omega_                 = Vec3::Zero();

    // Apply to all compiled structural elements. InertialLoad additionally
    // traverses ModelData::point_elements when this flag is enabled.
    inertia_load.region_                = model_data.elem_sets.all();
    inertia_load.consider_point_masses_ = consider_point_masses;

    // The InertialLoad implementation owns the established FEMaster sign
    // convention. The resulting field therefore already contains the signed
    // balancing inertia forces.
    inertia_load.apply(model_data, inertial_loads, Precision(0));

    // Superimpose the balancing inertia loads onto the externally applied loads
    global_load_mat += inertial_loads;

    // ---------------------------------------------------------------------
    // Verify the resulting rigid-body equilibrium
    // ---------------------------------------------------------------------

    const auto [resulting_force, resulting_moment] =
        compute_resultants(global_load_mat, center_of_gravity);

    logging::error(resulting_force.allFinite() && resulting_moment.allFinite(),
        "InertiaRelief: resulting force or moment contains NaN or Inf");

    logging::info(true, "resulting F/M:", std::setprecision(log_precision),
        std::setw(log_width), resulting_force(0),
        std::setw(log_width), resulting_force(1),
        std::setw(log_width), resulting_force(2),
        std::setw(log_width), resulting_moment(0),
        std::setw(log_width), resulting_moment(1),
        std::setw(log_width), resulting_moment(2));

    // Scale the admissible final residuals with the corresponding original
    // resultant magnitudes:
    //
    //     tol_F = max(tol_abs, tol_rel [1 + max |F_i|]),
    //     tol_M = max(tol_abs, tol_rel [1 + max |M_i|]).
    //
    // Force and moment use separate scales because they have different units
    // and may differ significantly in magnitude.
    const Precision force_tolerance = std::max(
        absolute_balance_tolerance,
        relative_balance_tolerance *
        (Precision(1) + external_force.cwiseAbs().maxCoeff())
    );

    const Precision moment_tolerance = std::max(
        absolute_balance_tolerance,
        relative_balance_tolerance *
        (Precision(1) + external_moment.cwiseAbs().maxCoeff())
    );

    // The infinity norm checks every Cartesian rigid-body component
    // independently and prevents one large residual from being hidden by the
    // remaining components.
    logging::error(resulting_force. lpNorm<Eigen::Infinity>() <= force_tolerance
                && resulting_moment.lpNorm<Eigen::Infinity>() <= moment_tolerance,
        "InertiaRelief: residual balance too large (|F|=", resulting_force.norm(),
        ", |M|=", resulting_moment.norm(), ")");
}

} // namespace fem
