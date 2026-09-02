// -----------------------------------------------------------------------------
// Hardening law
// -----------------------------------------------------------------------------

/**
 * @brief Builds the stress-scale tolerance used by constitutive local solves.
 *
 * The tolerance is intentionally scaled with a characteristic stress rather than
 * using a fixed absolute value. This keeps the local convergence criterion useful
 * for materials expressed in different unit systems.
 *
 *     tol = 1e-10 max(1, |scale_a|, |scale_b|)
 *
 * @param scale_a First characteristic stress.
 * @param scale_b Second characteristic stress.
 * @return Absolute stress tolerance.
 */
Precision stress_tolerance(Precision scale_a, Precision scale_b) {
    return Precision(1e-10)
         * std::max({Precision(1), std::abs(scale_a), std::abs(scale_b)});
}

/**
 * @brief Returns the initial yield stress of the tabulated hardening law.
 *
 * The first hardening point is required to correspond to zero equivalent plastic
 * strain and therefore defines the initial elastic domain.
 *
 * @param yield_curve Piecewise-linear isotropic hardening table.
 * @return Initial yield stress.
 */
Precision initial_yield_stress(const YieldCurve& yield_curve) {
    logging::error(!yield_curve.empty(),
        "J2: at least one yield point must be added before evaluation");

    return yield_curve.front().yield_stress;
}

/**
 * @brief Evaluates the piecewise-linear isotropic hardening law.
 *
 * For two adjacent tabulated points
 *
 *     (alpha_0, sigma_0), (alpha_1, sigma_1)
 *
 * the interpolation parameter and current yield stress are
 *
 *     xi      = (alpha - alpha_0) / (alpha_1 - alpha_0)
 *     sigma_y = sigma_0 + xi (sigma_1 - sigma_0).
 *
 * Values below the first point use the initial yield stress. Values beyond the
 * final point use the final tabulated stress, i.e. the continuation is perfectly
 * plastic.
 *
 * @param yield_curve Piecewise-linear isotropic hardening table.
 * @param eqp Accumulated equivalent plastic strain alpha.
 * @return Current yield stress sigma_y(alpha).
 */
Precision yield_stress_at(const YieldCurve& yield_curve, Precision eqp) {
    logging::error(!yield_curve.empty(),
        "J2: yield curve must not be empty");

    // The first entry also covers alpha = 0 and any harmless round-off below zero.
    if (yield_curve.size() == 1
        || eqp <= yield_curve.front().equivalent_plastic_strain) {
        return yield_curve.front().yield_stress;
    }

    // Search the active linear interpolation interval.
    for (std::size_t i = 1; i < yield_curve.size(); ++i) {
        const YieldPoint& left  = yield_curve[i - 1];
        const YieldPoint& right = yield_curve[i];

        if (eqp <= right.equivalent_plastic_strain) {
            const Precision strain_span = right.equivalent_plastic_strain
                                        - left.equivalent_plastic_strain;
            const Precision xi = (eqp - left.equivalent_plastic_strain) / strain_span;

            return left.yield_stress
                 + xi * (right.yield_stress - left.yield_stress);
        }
    }

    // The hardening curve is constant after its final tabulated point.
    return yield_curve.back().yield_stress;
}

/**
 * @brief Returns the active slope of the piecewise-linear hardening law.
 *
 * Inside one table interval the derivative is constant,
 *
 *     H = d sigma_y / d alpha
 *       = (sigma_1 - sigma_0) / (alpha_1 - alpha_0).
 *
 * At an interior table knot the hardening law is not differentiable. The interval
 * to the right is selected deliberately so that the local Newton system and the
 * algorithmic tangent use the same one-sided constitutive branch. Beyond the last
 * tabulated point H = 0 because the continuation is perfectly plastic.
 *
 * @param yield_curve Piecewise-linear isotropic hardening table.
 * @param eqp Accumulated equivalent plastic strain alpha.
 * @return Active hardening modulus H.
 */
Precision hardening_slope_at(const YieldCurve& yield_curve, Precision eqp) {
    if (yield_curve.size() < 2
        || eqp >= yield_curve.back().equivalent_plastic_strain) {
        return Precision(0);
    }

    for (std::size_t i = 1; i < yield_curve.size(); ++i) {
        const YieldPoint& left  = yield_curve[i - 1];
        const YieldPoint& right = yield_curve[i];

        if (eqp < right.equivalent_plastic_strain) {
            return (right.yield_stress - left.yield_stress)
                 / (right.equivalent_plastic_strain - left.equivalent_plastic_strain);
        }
    }

    return Precision(0);
}

/**
 * @brief Solves the scalar consistency equation of the small-strain radial return.
 *
 * For associative J2 plasticity with isotropic hardening the plastic multiplier
 * increment satisfies
 *
 *     r(dgamma) = q_trial
 *               - 3 G dgamma
 *               - sigma_y(alpha_n + dgamma)
 *               = 0.
 *
 * The hardening law is monotone but only piecewise differentiable. A bracketed
 * bisection is therefore deliberately used here: it is inexpensive for one scalar
 * unknown, cannot jump across a table kink, and is completely robust for perfect
 * plasticity as well as linear hardening.
 *
 * @param q_trial Trial von-Mises equivalent stress.
 * @param alpha_old Committed accumulated equivalent plastic strain.
 * @param shear Elastic shear modulus G.
 * @param yield_curve Piecewise-linear isotropic hardening table.
 * @return Plastic multiplier increment dgamma.
 */
Precision solve_small_plastic_increment(Precision q_trial,
                                        Precision alpha_old,
                                        Precision shear,
                                        const YieldCurve& yield_curve) {
    const Precision f_trial = q_trial - yield_stress_at(yield_curve, alpha_old);

    // The caller normally enters this function only for a plastic trial state.
    // Keeping the elastic guard makes the scalar solve safe when called directly.
    if (f_trial <= Precision(0)) {
        return Precision(0);
    }

    // The lower bound is the elastic trial state. Ignoring hardening gives the
    // classical perfect-plastic radial-return increment and therefore a valid
    // upper bound for a non-decreasing hardening curve.
    Precision lower = Precision(0);
    Precision upper = f_trial / (Precision(3) * shear);

    const Precision tolerance = stress_tolerance(
        Precision(3) * shear,
        initial_yield_stress(yield_curve)
    );

    for (Index iteration = 0; iteration < 80; ++iteration) {
        const Precision dgamma = Precision(0.5) * (lower + upper);
        const Precision residual = q_trial
                                 - Precision(3) * shear * dgamma
                                 - yield_stress_at(yield_curve, alpha_old + dgamma);

        if (std::abs(residual) <= tolerance) {
            return dgamma;
        }

        // Positive residual means the current plastic correction is too small;
        // negative residual means it is too large.
        if (residual > Precision(0)) {
            lower = dgamma;
        } else {
            upper = dgamma;
        }
    }

    // Eighty bisections are far beyond double-precision requirements for this
    // scalar problem. Returning the final bracket midpoint is therefore safer
    // than introducing a second failure mode solely for round-off saturation.
    return Precision(0.5) * (lower + upper);
}
