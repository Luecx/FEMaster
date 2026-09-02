// -----------------------------------------------------------------------------
// Small-strain return map and consistent tangent
// -----------------------------------------------------------------------------

/**
 * @brief Complete result of one small-strain J2 return map.
 *
 * Besides the converged Cauchy stress, the structure retains exactly the trial
 * quantities required by the analytical algorithmic tangent. None of these
 * values are persistent material history; the physical history is written only
 * into the seven-component material state vector.
 */
struct SmallResponse {
    VolumeStress stress;
    VolumeStress deviatoric_trial_stress;
    Mat3 flow_direction = Mat3::Zero();
    Precision q_trial   = Precision(0);
    Precision alpha_old = Precision(0);
    Precision dgamma    = Precision(0);
    bool plastic        = false;
};

/**
 * @brief Integrates the infinitesimal associative J2 model from one committed state.
 *
 * The routine uses the classical radial-return construction. The committed plastic
 * metric Cp is first converted to the logarithmic plastic strain used by the
 * infinitesimal model,
 *
 *     ep_n = 1/2 log(Cp_n).
 *
 * From the supplied total strain e the elastic trial strain and trial stress are
 *
 *     ee_trial    = e - ep_n
 *     sigma_trial = K tr(ee_trial) I + 2 G dev(ee_trial).
 *
 * The J2 trial measure is
 *
 *     s_trial = dev(sigma_trial)
 *     q_trial = sqrt(3/2 s_trial : s_trial)
 *     f_trial = q_trial - sigma_y(alpha_n).
 *
 * Elastic trial states are returned directly. Plastic trial states use a radial
 * correction with flow direction
 *
 *     N = 3/2 s_trial / q_trial,
 *
 * and the scalar consistency equation is solved for dgamma. The converged history
 * is then
 *
 *     ep_{n+1}    = ep_n + dgamma N
 *     Cp_{n+1}    = exp(2 ep_{n+1})
 *     alpha_{n+1} = alpha_n + dgamma.
 *
 * @param strain Generic six-component volume strain. The J2 kernel needs only the
 *               tensor/Voigt representation and is therefore independent of the
 *               public strain-measure wrapper.
 * @param state Working copy of the committed seven-component material state. The
 *              function overwrites it with the converged candidate state.
 * @param shear Elastic shear modulus G.
 * @param bulk Elastic bulk modulus K.
 * @param yield_curve Piecewise-linear isotropic hardening law.
 * @return Converged stress and tangent-relevant trial data.
 */
SmallResponse integrate_small_strain(const VolumeStrain& strain,
                                     State& state,
                                     Precision shear,
                                     Precision bulk,
                                     const YieldCurve& yield_curve) {
    SmallResponse result;

    // -------------------------------------------------------------------------
    // Reconstruct Cp_n explicitly from the persistent material state.
    //
    // State ordering:
    //     [Cp11, Cp22, Cp33, Cp23, Cp13, Cp12, alpha]
    //
    // Cp is symmetric, so only six independent entries are stored.
    // -------------------------------------------------------------------------
    Mat3 Cp_old;
    Cp_old << state[0], state[5], state[4],
              state[5], state[1], state[3],
              state[4], state[3], state[2];
    Cp_old = Precision(0.5) * (Cp_old + Cp_old.transpose());

    result.alpha_old = state[eqp_index];

    // -------------------------------------------------------------------------
    // Convert the plastic metric into the infinitesimal plastic strain.
    //
    //     Cp_n = exp(2 ep_n)
    //
    // therefore
    //
    //     ep_n = 1/2 log(Cp_n).
    //
    // Cp_n is symmetric positive definite. Its logarithm is therefore evaluated
    // spectrally from Cp_n = Q diag(lambda_i) Q^T.
    // -------------------------------------------------------------------------
    Eigen::SelfAdjointEigenSolver<Mat3> cp_solver(Cp_old);
    logging::error(cp_solver.info() == Eigen::Success,
        "J2: eigendecomposition failed for committed plastic metric");

    const Vec3 cp_eigenvalues = cp_solver.eigenvalues();
    const Precision cp_scale = std::max(
        Precision(1), cp_eigenvalues.cwiseAbs().maxCoeff()
    );
    logging::error(cp_eigenvalues.minCoeff()
                   > Precision(100) * std::numeric_limits<Precision>::epsilon() * cp_scale,
        "J2: committed plastic metric is not positive definite");

    Vec3 logarithmic_eigenvalues;
    logarithmic_eigenvalues(0) = std::log(cp_eigenvalues(0));
    logarithmic_eigenvalues(1) = std::log(cp_eigenvalues(1));
    logarithmic_eigenvalues(2) = std::log(cp_eigenvalues(2));

    Mat3 plastic_strain_old = Precision(0.5)
        * cp_solver.eigenvectors()
        * logarithmic_eigenvalues.asDiagonal()
        * cp_solver.eigenvectors().transpose();
    plastic_strain_old = Precision(0.5)
        * (plastic_strain_old + plastic_strain_old.transpose());

    // -------------------------------------------------------------------------
    // Build the elastic trial state.
    //
    // The VolumeStrain base class already owns the engineering-Voigt <-> tensor
    // conversion. Keeping that conversion in one place avoids duplicating shear
    // conventions in the constitutive model.
    // -------------------------------------------------------------------------
    Mat3 elastic_trial_strain = strain.tensor() - plastic_strain_old;
    elastic_trial_strain = Precision(0.5)
        * (elastic_trial_strain + elastic_trial_strain.transpose());

    const Precision elastic_trace = elastic_trial_strain.trace();
    const Mat3 elastic_deviator = elastic_trial_strain
        - (elastic_trace / Precision(3)) * Mat3::Identity();

    Mat3 trial_stress_tensor = bulk * elastic_trace * Mat3::Identity()
                             + Precision(2) * shear * elastic_deviator;
    trial_stress_tensor = Precision(0.5)
        * (trial_stress_tensor + trial_stress_tensor.transpose());

    // -------------------------------------------------------------------------
    // Evaluate the J2 yield function.
    //
    //     s_trial = sigma_trial - 1/3 tr(sigma_trial) I
    //     q_trial = sqrt(3/2 s_trial : s_trial)
    //
    // For matrices A and B the double contraction is simply the element-wise
    // product summed over all entries. It is intentionally written inline here
    // because the operation is used only as part of the J2 definition.
    // -------------------------------------------------------------------------
    const Mat3 deviatoric_trial_tensor = trial_stress_tensor
        - (trial_stress_tensor.trace() / Precision(3)) * Mat3::Identity();

    result.deviatoric_trial_stress = VolumeStress(deviatoric_trial_tensor);
    result.q_trial = std::sqrt(std::max(
        Precision(0),
        Precision(1.5)
            * (deviatoric_trial_tensor.array() * deviatoric_trial_tensor.array()).sum()
    ));

    const Precision yield_old = yield_stress_at(yield_curve, result.alpha_old);
    const Precision f_trial   = result.q_trial - yield_old;

    result.plastic = f_trial
        > stress_tolerance(Precision(3) * bulk, initial_yield_stress(yield_curve));

    // -------------------------------------------------------------------------
    // Elastic step.
    //
    // No history variable changes. The working state therefore already contains
    // the correct converged state and only the trial stress must be returned.
    // -------------------------------------------------------------------------
    if (!result.plastic) {
        result.stress = VolumeStress(trial_stress_tensor);
        return result;
    }

    logging::error(result.q_trial > Precision(0),
        "J2: positive yield residual with zero trial equivalent stress");

    // -------------------------------------------------------------------------
    // Plastic radial return.
    //
    // Associated J2 flow gives
    //
    //     N = df/dsigma
    //       = 3/2 s_trial / q_trial.
    //
    // The scalar consistency equation then supplies dgamma. Because the return is
    // radial, the direction N remains the trial direction throughout this update.
    // -------------------------------------------------------------------------
    result.flow_direction = (Precision(1.5) / result.q_trial)
                          * deviatoric_trial_tensor;

    result.dgamma = solve_small_plastic_increment(
        result.q_trial,
        result.alpha_old,
        shear,
        yield_curve
    );

    // -------------------------------------------------------------------------
    // Update the plastic strain and convert it back to the persistent metric Cp.
    //
    //     ep_{n+1} = ep_n + dgamma N
    //     Cp_{n+1} = exp(2 ep_{n+1}).
    //
    // 2 ep_{n+1} is symmetric, so its exponential is evaluated spectrally.
    // -------------------------------------------------------------------------
    Mat3 plastic_strain_new = plastic_strain_old
                            + result.dgamma * result.flow_direction;
    plastic_strain_new = Precision(0.5)
        * (plastic_strain_new + plastic_strain_new.transpose());

    const Mat3 exponent_argument = Precision(2) * plastic_strain_new;
    Eigen::SelfAdjointEigenSolver<Mat3> exp_solver(exponent_argument);
    logging::error(exp_solver.info() == Eigen::Success,
        "J2: eigendecomposition failed during small-strain plastic exponential");

    Vec3 exponential_eigenvalues;
    exponential_eigenvalues(0) = std::exp(exp_solver.eigenvalues()(0));
    exponential_eigenvalues(1) = std::exp(exp_solver.eigenvalues()(1));
    exponential_eigenvalues(2) = std::exp(exp_solver.eigenvalues()(2));

    Mat3 Cp_new = exp_solver.eigenvectors()
                * exponential_eigenvalues.asDiagonal()
                * exp_solver.eigenvectors().transpose();
    Cp_new = Precision(0.5) * (Cp_new + Cp_new.transpose());

    // Associated J2 plastic flow is isochoric, so det(Cp) = 1 analytically. The
    // normalization below removes only accumulated floating-point drift.
    const Precision determinant = Cp_new.determinant();
    logging::error(std::isfinite(determinant) && determinant > Precision(0),
        "J2: plastic metric lost positive determinant");
    Cp_new /= std::cbrt(determinant);
    Cp_new = Precision(0.5) * (Cp_new + Cp_new.transpose());

    // Store all seven physical state variables explicitly. This makes the state
    // layout visible at the update site and avoids hiding it behind generic copy
    // or packing helpers.
    state[0] = Cp_new(0, 0);
    state[1] = Cp_new(1, 1);
    state[2] = Cp_new(2, 2);
    state[3] = Cp_new(1, 2);
    state[4] = Cp_new(0, 2);
    state[5] = Cp_new(0, 1);
    state[6] = result.alpha_old + result.dgamma;

    // -------------------------------------------------------------------------
    // Correct the stress along the radial-return direction.
    //
    //     sigma_{n+1} = sigma_trial - 2 G dgamma N.
    // -------------------------------------------------------------------------
    Mat3 stress_tensor = trial_stress_tensor
                       - Precision(2) * shear * result.dgamma * result.flow_direction;
    stress_tensor = Precision(0.5) * (stress_tensor + stress_tensor.transpose());
    result.stress = VolumeStress(stress_tensor);

    return result;
}

/**
 * @brief Builds the analytically consistent small-strain J2 tangent.
 *
 * The returned matrix maps engineering-Voigt strain increments to physical-Voigt
 * Cauchy-stress increments. No constitutive perturbation is performed.
 *
 * In the elastic regime the tangent is simply isotropic elasticity. In the
 * plastic regime the scalar consistency equation gives
 *
 *     dgamma = dq_trial / (3 G + H),
 *
 * with active hardening modulus H. The normalized flow direction
 *
 *     N = 3/2 s_trial / q_trial
 *
 * is differentiated analytically as
 *
 *     dN = 3/(2 q_trial)
 *          [ds_trial - (dq_trial / q_trial) s_trial].
 *
 * The stress derivative then follows directly from
 *
 *     sigma = sigma_trial - 2 G dgamma N.
 *
 * @param base Converged return-map data from integrate_small_strain().
 * @param shear Elastic shear modulus G.
 * @param bulk Elastic bulk modulus K.
 * @param yield_curve Piecewise-linear isotropic hardening law.
 * @return Consistent six-by-six material tangent.
 */
Mat6 tangent_small(const SmallResponse& base,
                   Precision shear,
                   Precision bulk,
                   const YieldCurve& yield_curve) {
    Mat6 tangent = Mat6::Zero();

    // The hardening slope is needed only after plastic loading. At a table knot
    // hardening_slope_at() deliberately selects the right-sided active segment.
    const Precision hardening = base.plastic
        ? hardening_slope_at(yield_curve, base.alpha_old + base.dgamma)
        : Precision(0);

    const Mat3 deviatoric_trial_tensor = base.deviatoric_trial_stress.tensor();

    // Assemble the tangent one engineering-Voigt basis direction at a time. The
    // VolumeStrain class defines the engineering shear convention centrally, so
    // no manual factors of two are duplicated here.
    for (Index column = 0; column < 6; ++column) {
        VolumeStrain strain_direction;
        strain_direction.voigt().setZero();
        strain_direction.voigt()(column) = Precision(1);

        const Mat3 dstrain = strain_direction.tensor();
        const Precision dstrain_trace = dstrain.trace();
        const Mat3 dstrain_deviator = dstrain
            - (dstrain_trace / Precision(3)) * Mat3::Identity();

        // Elastic trial-stress derivative:
        //
        //     d sigma_trial = K tr(de) I + 2 G dev(de).
        Mat3 dstress = bulk * dstrain_trace * Mat3::Identity()
                     + Precision(2) * shear * dstrain_deviator;

        if (base.plastic) {
            // Only the deviatoric part enters q and the J2 flow direction.
            const Mat3 ddeviatoric_trial = dstress
                - (dstress.trace() / Precision(3)) * Mat3::Identity();

            // q = sqrt(3/2 s:s) gives dq = N:ds.
            const Precision dq =
                (base.flow_direction.array() * ddeviatoric_trial.array()).sum();

            // Differentiate the scalar consistency equation
            //
            //     q_trial - 3 G dgamma - sigma_y(alpha_n + dgamma) = 0
            //
            // to obtain
            //
            //     d(dgamma) = dq / (3 G + H).
            const Precision ddgamma = dq
                / (Precision(3) * shear + hardening);

            // Differentiate N = 3/2 s/q.
            const Mat3 dflow_direction = (Precision(1.5) / base.q_trial)
                * (ddeviatoric_trial
                   - (dq / base.q_trial) * deviatoric_trial_tensor);

            // Differentiate the corrected stress:
            //
            //     d sigma = d sigma_trial
            //             - 2 G [d(dgamma) N + dgamma dN].
            dstress -= Precision(2) * shear
                     * (ddgamma * base.flow_direction
                        + base.dgamma * dflow_direction);
        }

        // VolumeStress owns the physical-shear tensor -> Voigt conversion. The
        // column therefore has exactly the convention expected by FEMaster's
        // three-dimensional material tangent.
        tangent.col(column) = VolumeStress(dstress).voigt();
    }

    return tangent;
}
