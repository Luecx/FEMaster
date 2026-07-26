/**
 * @file rebalance_loads.cpp
 * @brief Implements rigid-body rebalancing of assembled nodal loads.
 *
 * The implementation removes the global resultant force and moment of an
 * assembled nodal load field without introducing support constraints or
 * computing rigid-body accelerations.
 *
 * The structural center of gravity is first determined from the model
 * elements. All applied nodal forces and, when available, nodal moments are
 * then reduced to a six-component resultant vector containing
 *
 *     (Fx, Fy, Fz, Mx, My, Mz).
 *
 * Six nodal force patterns are constructed:
 *
 * - uniform forces in the three global coordinate directions,
 * - rotational force patterns generated from the cross product between a
 *   global axis and the nodal position relative to the center of gravity.
 *
 * A dense 6x6 matrix maps the coefficients of these force patterns to their
 * global force and moment resultants. The coefficients required to cancel the
 * original resultants are obtained through a direct full-pivoting LU solve.
 * A damped least-squares solution is used when the resultant matrix is
 * singular or numerically unsuitable for the direct solve.
 *
 * The final correction is added only to the translational components of the
 * supplied load field. A postprocessing check verifies that the resulting
 * global force and moment satisfy the prescribed numerical tolerances.
 *
 * @see rebalance_loads
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#include "rebalance_loads.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"

#include <Eigen/Cholesky>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <utility>

namespace fem {

/**
 * Rebalances an assembled nodal load field in the six-dimensional rigid-body
 * load space.
 *
 * The routine first determines the structural center of gravity and evaluates
 * the global resultant force and moment of the supplied load field about this
 * point. Existing nodal moments are included when the field contains six or
 * more components.
 *
 * Six nodal force patterns are subsequently evaluated. Their global
 * resultants form the columns of a dense matrix `A`, such that
 *
 *     A * coefficients = applied_resultants.
 *
 * The coefficient vector is chosen to satisfy
 *
 *     A * coefficients = -initial_resultants.
 *
 * A full-pivoting LU decomposition is preferred because the system contains
 * only six equations. If the matrix is singular or the direct solution is not
 * finite, the routine solves the damped normal equations
 *
 *     (A^T A + lambda I) * coefficients = -A^T * initial_resultants.
 *
 * The resulting correction consists exclusively of nodal forces. No nodal
 * moments, constraints, displacements, accelerations or stiffness terms are
 * generated or modified.
 *
 * Degenerate geometries may not support all three independent rotational load
 * patterns. In that case, the damped least-squares fallback determines the
 * closest finite correction, and the final balance check verifies whether the
 * remaining residual is acceptable.
 *
 * @param model_data Global model data containing nodal positions and structural
 *                   elements used to determine the center of gravity.
 * @param global_load_mat Assembled nodal load field modified directly by the
 *                        routine. Components 0 to 2 contain nodal forces;
 *                        components 3 to 5 are interpreted as nodal moments
 *                        when available.
 */
void rebalance_loads(model::ModelData& model_data, model::Field& global_load_mat) {
    using Matrix6 = Eigen::Matrix<Precision, 6, 6>;
    using Vector6 = Eigen::Matrix<Precision, 6, 1>;

    // Numerical tolerances for the force and moment balance checks
    constexpr Precision relative_balance_tolerance = Precision(1e-10);
    constexpr Precision absolute_balance_tolerance = Precision(1e-12);

    // Regularisation parameters for the damped least-squares fallback
    constexpr Precision minimum_dls_regularisation  = Precision(1e-14);
    constexpr Precision relative_dls_regularisation = Precision(1e-12);

    // Formatting parameters for the diagnostic output
    constexpr int log_width     = 12;
    constexpr int log_precision = 3;

    // Validate the model fields required for the resultant and center-of-gravity
    // calculations
    logging::error(model_data.positions != nullptr,
        "RigidBodyRebalancing: positions not initialized");

    auto& positions = *model_data.positions;

    logging::error(global_load_mat.rows == positions.rows,
        "RigidBodyRebalancing: global_load_mat.rows != positions.rows (", global_load_mat.rows, " vs ", positions.rows, ")");
    logging::error(global_load_mat.components >= 3,
        "RigidBodyRebalancing: global_load_mat.components must be >= 3");

    /**
     * Computes the global force and moment resultants of a nodal load field.
     *
     * Nodal forces contribute directly to the force resultant and through their
     * lever arms to the moment resultant. If the load field contains at least
     * six components, the nodal moment components are added directly to the
     * global moment.
     *
     * All moments are evaluated about the supplied global reference point.
     */
    const auto compute_resultants = [&](const model::Field& loads, const Vec3& reference_point) -> std::pair<Vec3, Vec3> {
        Vec3 resultant_force  = Vec3::Zero();
        Vec3 resultant_moment = Vec3::Zero();

        const bool has_nodal_moments = loads.components >= 6;

        // Reduce every nodal load to the selected global reference point
        for (Index node = 0; node < loads.rows; ++node) {
            const Vec3 nodal_force{
                loads(node, 0),
                loads(node, 1),
                loads(node, 2)
            };

            // Nodal moments are optional because many load fields contain only
            // the three translational force components
            Vec3 nodal_moment = Vec3::Zero();

            if (has_nodal_moments) {
                nodal_moment = Vec3{
                    loads(node, 3),
                    loads(node, 4),
                    loads(node, 5)
                };
            }

            // Transform the nodal force into a moment about the selected
            // reference point
            const Vec3 position = positions.row_vec3(node);
            const Vec3 lever    = position - reference_point;

            resultant_force  += nodal_force;
            resultant_moment += nodal_moment + lever.cross(nodal_force);
        }

        return {resultant_force, resultant_moment};
    };

    // Determine the total structural mass and its first spatial moment. The
    // integration mode matches the convention used by the structural
    // inertia-relief implementation.
    Precision total_mass        = Precision(0);
    Vec3      first_mass_moment = Vec3::Zero();

    for (const auto& element : model_data.elements) {
        // Empty element entries do not contribute to the structural mass
        if (!element) {
            continue;
        }

        // Non-structural elements do not provide the volume or section
        // integration required for the center-of-gravity calculation
        auto* structural_element = element->as<model::StructuralElement>();

        if (!structural_element) {
            continue;
        }

        // Integrating one gives the element mass according to the selected
        // structural integration convention
        total_mass += structural_element->integrate_scalar_field(true,
            [](const Vec3&) { return Precision(1);});

        // Integrating the global position gives the first mass moment required
        // for the center of gravity
        first_mass_moment += structural_element->integrate_vector_field(true,
            [](const Vec3& position) { return position;});
    }

    // A finite positive structural mass is required to define the center of
    // gravity and therefore the reference point of the moment balance
    logging::error(total_mass > Precision(0) && std::isfinite(total_mass),
        "RigidBodyRebalancing: total mass is zero or invalid; cannot determine center of gravity");
    logging::error(first_mass_moment.allFinite(),
        "RigidBodyRebalancing: first mass moment contains NaN or Inf");

    const Vec3 center_of_gravity = first_mass_moment / total_mass;

    // Evaluate the initial force and moment imbalance about the structural
    // center of gravity
    const auto [initial_force, initial_moment] = compute_resultants(global_load_mat, center_of_gravity);

    logging::error(initial_force.allFinite() && initial_moment.allFinite(),
        "RigidBodyRebalancing: initial resultants contain NaN or Inf");

    // Scale the admissible residuals with the magnitude of the initial force
    // and moment while retaining an absolute lower tolerance
    const Precision force_tolerance = std::max(
        absolute_balance_tolerance,
        relative_balance_tolerance * (Precision(1) + initial_force.cwiseAbs().maxCoeff())
    );

    const Precision moment_tolerance = std::max(
        absolute_balance_tolerance,
        relative_balance_tolerance * (Precision(1) + initial_moment.cwiseAbs().maxCoeff())
    );

    // Print the initial six-component rigid-body load resultant
    logging::info(true, "RigidBodyRebalancing:");
    logging::info(true, "              ",
        std::setw(log_width), "x",
        std::setw(log_width), "y",
        std::setw(log_width), "z",
        std::setw(log_width), "rx",
        std::setw(log_width), "ry",
        std::setw(log_width), "rz");
    logging::info(true, "initial F/M  :", std::setprecision(log_precision),
        std::setw(log_width), initial_force(0),
        std::setw(log_width), initial_force(1),
        std::setw(log_width), initial_force(2),
        std::setw(log_width), initial_moment(0),
        std::setw(log_width), initial_moment(1),
        std::setw(log_width), initial_moment(2));

    // Avoid constructing and solving the correction system when the load field
    // already satisfies the required rigid-body equilibrium
    if (initial_force.lpNorm<Eigen::Infinity>() <= force_tolerance &&
        initial_moment.lpNorm<Eigen::Infinity>() <= moment_tolerance) {
        logging::info(true, "RigidBodyRebalancing: already balanced within tolerance.");
        return;
    }

    // Uniform nodal weights distribute each correction pattern equally among
    // all model nodes.
    //
    // TODO: Replace the uniform weights with lumped nodal masses if the
    // correction should follow the physical mass distribution.
    const DynamicVector basis_weights = DynamicVector::Ones(positions.rows);

    /**
     * Evaluates one of the six rigid-load basis forces at a selected node.
     *
     * Basis indices 0 to 2 define uniform nodal forces in the global x, y and z
     * directions. Basis indices 3 to 5 define rotational force patterns
     *
     *     force_i = weight_i * (axis x radius_i),
     *
     * where `radius_i` is measured from the structural center of gravity.
     *
     * A rotational pattern is not assumed to produce an exactly pure moment.
     * Any translational or rotational coupling caused by the nodal geometry is
     * represented explicitly in the subsequently assembled 6x6 resultant
     * matrix.
     */
    const auto evaluate_basis_force = [&](Index node, Index basis) -> Vec3 {
        // Select the global axis associated with the current translation or
        // rotation basis
        Vec3 axis = Vec3::Zero();
        axis(basis % 3) = Precision(1);

        // Translational basis patterns apply an equally directed force to every
        // node
        if (basis < 3) {
            return basis_weights(node) * axis;
        }

        // Rotational basis patterns apply forces tangential to the selected
        // global rotation axis
        const Vec3 position = positions.row_vec3(node);
        const Vec3 radius   = position - center_of_gravity;

        return basis_weights(node) * axis.cross(radius);
    };

    /**
     * Computes the force and moment resultants of one rigid-load basis pattern.
     *
     * The basis is evaluated directly at every node so no temporary load field
     * needs to be allocated.
     */
    const auto compute_basis_resultants = [&](Index basis) -> std::pair<Vec3, Vec3> {
        Vec3 resultant_force  = Vec3::Zero();
        Vec3 resultant_moment = Vec3::Zero();

        // Reduce the complete nodal basis pattern to the center of gravity
        for (Index node = 0; node < positions.rows; ++node) {
            const Vec3 position    = positions.row_vec3(node);
            const Vec3 radius      = position - center_of_gravity;
            const Vec3 nodal_force = evaluate_basis_force(node, basis);

            resultant_force  += nodal_force;
            resultant_moment += radius.cross(nodal_force);
        }

        return {resultant_force, resultant_moment};
    };

    // Assemble the matrix mapping the six basis coefficients to the global
    // force and moment resultants
    Matrix6 resultant_matrix = Matrix6::Zero();

    for (Index basis = 0; basis < 6; ++basis) {
        const auto [basis_force, basis_moment] = compute_basis_resultants(basis);

        resultant_matrix.block<3, 1>(0, basis) = basis_force;
        resultant_matrix.block<3, 1>(3, basis) = basis_moment;
    }

    // Collect the initial force and moment imbalance in a single right-hand
    // side vector
    Vector6 initial_resultants = Vector6::Zero();

    initial_resultants.segment<3>(0) = initial_force;
    initial_resultants.segment<3>(3) = initial_moment;

    // Solve for coefficients whose induced resultants cancel the initial
    // imbalance
    Vector6 coefficients    = Vector6::Zero();
    bool    direct_solution = false;

    {
        // Full-pivoting LU handles the small dense system and provides an
        // explicit invertibility check before solving
        Eigen::FullPivLU<Matrix6> decomposition(resultant_matrix);

        if (decomposition.isInvertible()) {
            coefficients    = decomposition.solve(-initial_resultants);
            direct_solution = coefficients.allFinite();
        }
    }

    // Use damped least squares when the six load patterns do not form an
    // invertible resultant basis for the current nodal geometry
    if (!direct_solution) {
        // Build the symmetric normal-equation matrix
        const Matrix6 normal_matrix = resultant_matrix.transpose() * resultant_matrix;

        // Scale the damping with the largest diagonal coefficient while
        // retaining a strictly positive absolute lower bound
        const Precision diagonal_scale = normal_matrix.diagonal().cwiseAbs().maxCoeff();

        const Precision regularisation = std::max(
            minimum_dls_regularisation,
            relative_dls_regularisation * (Precision(1) + diagonal_scale)
        );

        // Add the regularisation directly to the diagonal instead of
        // constructing a separate identity matrix
        Matrix6 regularised_matrix = normal_matrix;
        regularised_matrix.diagonal().array() += regularisation;

        const Vector6 right_hand_side = -resultant_matrix.transpose() * initial_resultants;

        // The regularised normal matrix is symmetric and is solved through an
        // LDLT factorisation
        Eigen::LDLT<Matrix6> decomposition(regularised_matrix);

        logging::error(decomposition.info() == Eigen::Success,
            "RigidBodyRebalancing: damped least-squares LDLT factorization failed");

        coefficients = decomposition.solve(right_hand_side);

        logging::error(decomposition.info() == Eigen::Success,
            "RigidBodyRebalancing: damped least-squares solve failed");
        logging::error(coefficients.allFinite(),
            "RigidBodyRebalancing: damped least-squares coefficients contain NaN or Inf");
    }

    // Evaluate the force and moment resultants produced by the selected
    // coefficients for diagnostic output
    const Vector6 applied_resultants = resultant_matrix * coefficients;
    const Vec3    applied_force      = applied_resultants.segment<3>(0);
    const Vec3    applied_moment     = applied_resultants.segment<3>(3);

    logging::info(true, "coeff (Tx,Ty,Tz,Rx,Ry,Rz):",
        std::setprecision(log_precision),
        std::setw(log_width), coefficients(0),
        std::setw(log_width), coefficients(1),
        std::setw(log_width), coefficients(2),
        std::setw(log_width), coefficients(3),
        std::setw(log_width), coefficients(4),
        std::setw(log_width), coefficients(5));
    logging::info(true, "applied dF/M :",
        std::setprecision(log_precision),
        std::setw(log_width), applied_force(0),
        std::setw(log_width), applied_force(1),
        std::setw(log_width), applied_force(2),
        std::setw(log_width), applied_moment(0),
        std::setw(log_width), applied_moment(1),
        std::setw(log_width), applied_moment(2));

    // Construct and apply the nodal force correction. The basis definition is
    // reused here so the applied load field remains exactly consistent with
    // the matrix used to determine the coefficients.
    for (Index node = 0; node < positions.rows; ++node) {
        Vec3 correction_force = Vec3::Zero();

        // Superimpose all six rigid-load basis patterns according to the solved
        // coefficient vector
        for (Index basis = 0; basis < 6; ++basis) {
            correction_force += coefficients(basis) * evaluate_basis_force(node, basis);
        }

        // Modify only the translational nodal load components
        global_load_mat(node, 0) += correction_force(0);
        global_load_mat(node, 1) += correction_force(1);
        global_load_mat(node, 2) += correction_force(2);
    }

    // Recompute the actual resultants from the modified load field rather than
    // relying only on the small dense correction system
    const auto [resulting_force, resulting_moment] = compute_resultants(global_load_mat, center_of_gravity);

    logging::info(true, "resulting F/M:",
        std::setprecision(log_precision),
        std::setw(log_width), resulting_force(0),
        std::setw(log_width), resulting_force(1),
        std::setw(log_width), resulting_force(2),
        std::setw(log_width), resulting_moment(0),
        std::setw(log_width), resulting_moment(1),
        std::setw(log_width), resulting_moment(2));

    // Require both the force and moment residuals to satisfy their scaled
    // infinity-norm tolerances
    logging::error(resulting_force. lpNorm<Eigen::Infinity>() <= force_tolerance
                && resulting_moment.lpNorm<Eigen::Infinity>() <= moment_tolerance,
        "RigidBodyRebalancing: residual balance too large (|F|=", resulting_force.norm(),
        ", |M|=", resulting_moment.norm(), ")");
}

} // namespace fem