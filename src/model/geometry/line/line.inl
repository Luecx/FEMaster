/**
 * @file line.inl
 * @brief Implements common geometry operations for templated line elements.
 *
 * The routines gather nodal coordinates, evaluate physical mappings, integrate
 * line length and solve the closest-point projection problem. Element-specific
 * shape functions and derivatives remain implemented by the concrete line
 * classes.
 *
 * @see Line
 * @see LineInterface
 */

namespace fem {
namespace model {

template<Index N, IsoParametricLineRange CR>
Line<N, CR>::Line(const std::array<ID, N>& node_ids)
    : LineInterface(N),
      node_ids    (node_ids) {}

template<Index N, IsoParametricLineRange CR>
const ID* Line<N, CR>::nodes() const {
    return node_ids.data();
}

template<Index N, IsoParametricLineRange CR>
ID* Line<N, CR>::nodes() {
    return node_ids.data();
}

template<Index N, IsoParametricLineRange CR>
DynamicVector Line<N, CR>::shape_function(Precision r) const {
    // Evaluate the fixed-size element interpolation through the common dynamic
    // interface used by generic line geometry code.
    const StaticMatrix<N, 1> shape_func = this->_shape_function(r);
    return DynamicVector{shape_func};
}

/**
 * Gathers the global coordinates of all nodes attached to the line.
 *
 * @param node_coords_system Global nodal coordinate field.
 * @return Matrix containing one global coordinate row per line node.
 */
template<Index N, IsoParametricLineRange CR>
StaticMatrix<N, 3> Line<N, CR>::node_coords_global(const Field& node_coords_system) const {
    StaticMatrix<N, 3> res {};

    // Gather the global coordinate row belonging to every local line node.
    for (Index i = 0; i < N; i++) {
        const Index row = static_cast<Index>(node_ids[i]);
        res.row(i) = node_coords_system.row_vec3(row).transpose();
    }

    return res;
}

/**
 * Computes the physical length of the line by numerical integration.
 *
 * The physical tangent is assembled from the nodal coordinates and shape
 * function derivatives. Its norm is integrated over the natural interval
 * selected by `CR`.
 *
 * @param node_coords_system Global nodal coordinate field.
 * @return Physical length of the line element.
 */
template<Index N, IsoParametricLineRange CR>
Precision Line<N, CR>::length(const Field& node_coords_system) const {
    using namespace fem::math::quadrature;

    // Select the quadrature domain matching the natural-coordinate range.
    Domain domain = (CR == MINUS_ONE_TO_ONE) ? DOMAIN_ISO_LINE_A : DOMAIN_ISO_LINE_B;

    // Use a quadratic-order rule to integrate the physical tangent norm.
    Quadrature quadrature(domain, ORDER_QUADRATIC);

    // Gather nodal coordinates once because they are reused at every point.
    const StaticMatrix<N, 3> node_coords_global = this->node_coords_global(node_coords_system);

    // Integrate the physical tangent norm over the complete natural domain.
    return quadrature.integrate([&](Precision r, Precision, Precision) -> Precision {
        Vec3 dx_dr = Vec3::Zero();
        const StaticMatrix<N, 1> dN_dr = this->shape_derivative(r);

        // Assemble the physical tangent from nodal coordinates and shape
        // function derivatives.
        for (Index i = 0; i < N; i++) {
            dx_dr += dN_dr(i) * node_coords_global.row(i);
        }

        return dx_dr.norm();
    });
}

/**
 * Finds the natural coordinate of the closest point on the line.
 *
 * Newton iterations minimize one half of the squared distance between the
 * supplied global point and the isoparametric line. Several topology-dependent
 * initial guesses are evaluated, and the closest converged candidate is kept.
 * When `clip` is enabled, every iterate is limited to the natural interval.
 *
 * @param p Global point to project onto the line.
 * @param node_coords_system Global nodal coordinate field.
 * @param clip Whether Newton iterates must remain inside the natural range.
 * @return Natural coordinate of the closest projection found.
 */
template<Index N, IsoParametricLineRange CR>
Precision Line<N, CR>::global_to_local(
    const Vec3&  p,
    const Field& node_coords_system,
    bool         clip
) const {
    constexpr int       max_iter = 32;
    constexpr Precision eps      = 1e-12;
    std::vector<Precision> initial_guesses;

    // Use topology-dependent initial guesses to reduce the risk of converging
    // to an unsuitable stationary point on a curved line.
    if (N == 3) {
        initial_guesses = {min_r(), max_r()};
    } else if (N > 3) {
        initial_guesses = {min_r(), (min_r() + max_r()) / Precision(2.0), max_r()};
    } else {
        initial_guesses = {(min_r() + max_r()) / Precision(2.0)};
    }

    // Gather nodal coordinates once for all initial guesses and iterations.
    const auto node_coords_global = this->node_coords_global(node_coords_system);

    Precision best_r               = min_r();
    Precision min_distance_squared = std::numeric_limits<Precision>::max();

    // Minimize the squared-distance objective from every initial guess.
    for (const Precision initial_r : initial_guesses) {
        Precision r = initial_r;

        // Apply Newton's method to the one-dimensional projection problem.
        for (Index iter = 0; iter < max_iter; iter++) {
            const StaticMatrix<N, 1> N_vals  = _shape_function(r);
            const StaticMatrix<N, 1> dN_dr   = shape_derivative(r);
            const StaticMatrix<N, 1> d2N_dr2 = shape_second_derivative(r);

            Vec3 x_r      = Vec3::Zero();
            Vec3 dx_dr    = Vec3::Zero();
            Vec3 d2x_dr2 = Vec3::Zero();

            // Evaluate the physical position, tangent and curvature at the
            // current natural coordinate.
            for (Index i = 0; i < N; i++) {
                x_r      += N_vals(i)  * node_coords_global.row(i);
                dx_dr    += dN_dr(i)   * node_coords_global.row(i);
                d2x_dr2 += d2N_dr2(i) * node_coords_global.row(i);
            }

            // Evaluate the gradient and Hessian of one half of the squared
            // distance objective.
            const Vec3      diff    = x_r - p;
            const Precision dD_dr   = diff.dot(dx_dr);
            const Precision d2D_dr2 = dx_dr.dot(dx_dr) + diff.dot(d2x_dr2);

            // A nearly singular Hessian cannot provide a stable Newton step.
            if (std::abs(d2D_dr2) < 1e-14) {
                break;
            }

            // Apply the Newton correction and optionally enforce the natural
            // coordinate bounds.
            const Precision delta_r = -dD_dr / d2D_dr2;
            r += delta_r;

            if (clip) {
                if (r < min_r()) r = min_r();
                if (r > max_r()) r = max_r();
            }

            // Stop once the natural-coordinate correction is sufficiently small.
            if (std::abs(delta_r) < eps) {
                break;
            }
        }

        // Compare the converged candidate using squared physical distance.
        const Vec3      x_r              = this->local_to_global(r, node_coords_system);
        const Precision distance_squared = (x_r - p).squaredNorm();

        if (distance_squared < min_distance_squared) {
            min_distance_squared = distance_squared;
            best_r               = r;
        }
    }

    return best_r;
}

/**
 * Maps a natural line coordinate to its physical global position.
 *
 * @param r Natural line coordinate.
 * @param node_coords_system Global nodal coordinate field.
 * @return Interpolated global position.
 */
template<Index N, IsoParametricLineRange CR>
Vec3 Line<N, CR>::local_to_global(Precision r, const Field& node_coords_system) const {
    // Gather nodal coordinates and evaluate the shape functions at r.
    const auto node_coords_global = this->node_coords_global(node_coords_system);
    const auto N_vals             = _shape_function(r);
    Vec3       res                = Vec3::Zero();

    // Assemble the isoparametric position from all nodal contributions.
    for (Index i = 0; i < N; i++) {
        res += N_vals(i) * node_coords_global.row(i);
    }

    return res;
}

template<Index N, IsoParametricLineRange CR>
constexpr Precision Line<N, CR>::min_r() const {
    return CR == MINUS_ONE_TO_ONE ? -1 : 0;
}

template<Index N, IsoParametricLineRange CR>
constexpr Precision Line<N, CR>::max_r() const {
    return 1;
}

template<Index N, IsoParametricLineRange CR>
const std::array<ID, N>& Line<N, CR>::get_node_ids() const {
    return node_ids;
}

} // namespace model
} // namespace fem
