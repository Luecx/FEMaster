#include "optimise.h"

namespace fem {
namespace optimise {

/**
 * @brief Perform optimization to find the closest point in isoparametric coordinates (r, s).
 *
 * This function finds the (r, s) coordinates on an isoparametric element surface that are closest
 * to a given 3D point. Depending on the domain type, the optimization is either constrained
 * or unconstrained. The optimization uses the Adam optimizer to find the optimal parameters.
 *
 * @param point_3d Target 3D point for which the closest isoparametric coordinates are to be found.
 * @param surface_func Lambda function that takes r, s and returns a 3D point on the surface.
 * @param domain ElementDomain specifying the type of the domain (DOMAIN_ISO_TRI or DOMAIN_ISO_QUAD).
 * @param constrained Boolean flag to specify if constraints should be applied (default: true).
 * @param max_distance Maximum allowed distance for optimization. If the distance exceeds this value,
 *                     the function will terminate and return the initial guess.
 * @return StaticVector<2> Optimized r, s coordinates in the isoparametric domain.
 */
StaticVector<2> find_isoparametric_point(StaticVector<3>                                      point_3d,
                                         std::function<StaticVector<3>(Precision, Precision)> surface_func,
                                         fem::quadrature::Domain                              domain,
                                         bool                                                 constrained) {

    // Set Adam optimizer parameters
    Precision alpha     = 0.1;    // Step size (learning rate)
    Precision beta1     = 0.9;    // Exponential decay rate for the first moment
    Precision beta2     = 0.999;  // Exponential decay rate for the second moment
    Precision epsilon   = 1e-8;   // Small constant to avoid division by zero

    Precision tolerance = 1e-6;   // Tolerance for change in parameters
    Precision delta_t   = 1e-6;   // Small step size for numerical differentiation
    int max_iter        = 1000;   // Maximum number of iterations

    int iter = 0;

    // Initialize first and second moment vectors for Adam
    StaticVector<2> m = StaticVector<2>::Zero();
    StaticVector<2> v = StaticVector<2>::Zero();

    // Initial guess for r, s
    StaticVector<2> r_s;
    if (domain == fem::quadrature::Domain::DOMAIN_ISO_TRI) {
        r_s = StaticVector<2> {0.25, 0.25};
    } else {
        r_s = StaticVector<2> {0.0, 0.0};
    }

    for (iter = 0; iter < max_iter; ++iter) {
        // Get the surface point based on the current r, s
        StaticVector<3> surface_point = surface_func(r_s[0], r_s[1]);

        // Compute the distance and gradient using numerical differentiation
        StaticVector<3> diff     = surface_point - point_3d;
        Precision       distance = diff.norm();

        // Calculate numerical gradient
        StaticVector<2> gradient;
        for (int i = 0; i < 2; ++i) {
            StaticVector<2> r_s_plus = r_s;
            r_s_plus[i] += delta_t;
            StaticVector<3> surface_point_plus = surface_func(r_s_plus[0], r_s_plus[1]);
            gradient[i]                        = (surface_point_plus - surface_point).dot(diff) / delta_t;
        }

        // Update first and second moments
        m = beta1 * m + (1.0 - beta1) * gradient;
        v = beta2 * v + (1.0 - beta2) * gradient.cwiseProduct(gradient);

        // Use vector epsilon to match the type of v_hat
        StaticVector<2> epsilon_vec = StaticVector<2>::Constant(epsilon);

        // Compute Adam step
        StaticVector<2> delta_r_s = alpha * m.cwiseQuotient(v.cwiseSqrt() + epsilon_vec);

        // Update r, s using Adam step
        StaticVector<2> r_s_new = r_s - delta_r_s;

        // Track total change, including clipping adjustments
        StaticVector<2> r_s_change = delta_r_s;

        // Apply constraints if needed and track clipping changes
        if (constrained) {
            if (domain == fem::quadrature::Domain::DOMAIN_ISO_TRI) {
                // Triangle domain constraints: 0 <= r <= 1, 0 <= s <= 1, r + s <= 1
                r_s_new[0] = std::max(0.0, std::min(1.0, r_s_new[0]));
                r_s_new[1] = std::max(0.0, std::min(1.0, r_s_new[1]));
                if (r_s_new[0] + r_s_new[1] > 1.0) {
                    Precision diff = r_s_new[0] + r_s_new[1] - 1.0;
                    // Adjust both r_s_new[0] and r_s_new[1] to satisfy the constraint
                    r_s_new[0] -= diff / 2.0;
                    r_s_new[1] -= diff / 2.0;
                }
            } else if (domain == fem::quadrature::Domain::DOMAIN_ISO_QUAD) {
                // Quadrilateral domain constraints: -1 <= r <= 1, -1 <= s <= 1
                r_s_new[0] = std::max(-1.0, std::min(1.0, r_s_new[0]));
                r_s_new[1] = std::max(-1.0, std::min(1.0, r_s_new[1]));
            }
        }

        // Compute total change after applying constraints
        r_s_change = r_s_new - r_s;

        // Check for convergence based on the change in parameters including clipping effects
        if (r_s_change.norm() < tolerance) {
            std::cout << "Converged after " << iter << " iterations." << std::endl;
            break;
        }

        // Update r_s to the new value
        r_s = r_s_new;
    }

    std::cout << "Optimization terminated after " << iter << " iterations using Adam optimizer." << std::endl;

    return r_s;
}

}  // namespace optimise
}  // namespace fem
