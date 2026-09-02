// -----------------------------------------------------------------------------
// Shell and axial stress reductions
// -----------------------------------------------------------------------------

/**
 * @brief Solves the shell thickness strain required by the plane-stress condition.
 *
 * Integrated shell material laws receive five in-plane/transverse-shear strain
 * components. The missing thickness-normal strain e33 is determined from
 *
 *     S33(e33) = 0.
 *
 * A scalar Newton iteration is used. Crucially, the Newton derivative is not
 * approximated numerically: it is the (33,33) entry of the same consistent
 * three-dimensional J2 tangent returned by the current constitutive candidate,
 *
 *     de33 = -S33 / C3333.
 *
 * The complete three-dimensional strain vector uses FEMaster's standard volume
 * ordering
 *
 *     [xx, yy, zz, gamma_yz, gamma_xz, gamma_xy].
 *
 * The shell input ordering is mapped explicitly inside this function so the
 * reduction can be understood without consulting another helper.
 *
 * @tparam Finite false for linearized/Cauchy response, true for
 *                Green-Lagrange/PK2 response.
 * @param shell_strain Five shell material strain components.
 * @param committed Seven-component committed J2 state.
 * @param converged_state Receives the state associated with the converged e33.
 * @param youngs Young's modulus E, used only for scaling the local tolerance.
 * @param poisson Poisson ratio nu, used for the elastic initial guess.
 * @param shear Elastic shear modulus G.
 * @param bulk Elastic bulk modulus K.
 * @param yield_curve Piecewise-linear isotropic hardening law.
 * @param converged_tangent Receives the full three-dimensional consistent tangent
 *                          before plane-stress condensation.
 * @return Converged full three-dimensional stress.
 */
template<bool Finite>
VolumeStress solve_shell_plane_stress(const Vec5& shell_strain,
                                      const State& committed,
                                      State& converged_state,
                                      Precision youngs,
                                      Precision poisson,
                                      Precision shear,
                                      Precision bulk,
                                      const YieldCurve& yield_curve,
                                      Mat6& converged_tangent) {
    // -------------------------------------------------------------------------
    // Elastic plane-stress estimate for the missing thickness strain.
    //
    // For isotropic elasticity under sigma_33 = 0,
    //
    //     e33 = -nu / (1 - nu) (e11 + e22).
    //
    // This is only the Newton starting value; plasticity and finite kinematics are
    // handled by the complete three-dimensional material update below.
    // -------------------------------------------------------------------------
    Precision e33 = -poisson / (Precision(1) - poisson)
                  * (shell_strain(0) + shell_strain(1));

    if constexpr (Finite) {
        // Green-Lagrange strain must satisfy C = I + 2 E > 0. The lower clamp is
        // only a safe initial guess; every subsequent candidate is checked by its
        // full metric eigenvalues before it is accepted.
        e33 = std::max(e33, Precision(-0.49));
    }

    const Precision tolerance = stress_tolerance(
        youngs, initial_yield_stress(yield_curve)
    );

    VolumeStress stress;
    bool converged = false;

    for (Index iteration = 0; iteration < 30; ++iteration) {
        // ---------------------------------------------------------------------
        // Build the complete six-component strain state explicitly.
        //
        // Shell ordering used here:
        //     [e11, e22, gamma12, gamma13, gamma23]
        //
        // Volume ordering:
        //     [e11, e22, e33, gamma23, gamma13, gamma12].
        // ---------------------------------------------------------------------
        VolumeStrain volume_strain;
        volume_strain.voigt()(0) = shell_strain(0);
        volume_strain.voigt()(1) = shell_strain(1);
        volume_strain.voigt()(2) = e33;
        volume_strain.voigt()(3) = shell_strain(4);
        volume_strain.voigt()(4) = shell_strain(3);
        volume_strain.voigt()(5) = shell_strain(2);

        if constexpr (Finite) {
            // Finite-strain admissibility follows directly from
            //
            //     C = I + 2 E.
            //
            // A positive-definite C is required for a physical deformation.
            const Mat3 C = Mat3::Identity()
                         + Precision(2) * volume_strain.tensor();
            Eigen::SelfAdjointEigenSolver<Mat3> solver(
                Precision(0.5) * (C + C.transpose()),
                Eigen::EigenvaluesOnly
            );

            logging::error(solver.info() == Eigen::Success
                           && solver.eigenvalues().minCoeff() > Precision(1e-12),
                "J2: shell plane-stress iterate left admissible Green-Lagrange domain");
        }

        // ---------------------------------------------------------------------
        // Every plane-stress Newton candidate must start from exactly the same
        // committed material history. The state is therefore copied explicitly
        // component-by-component instead of relying on a hidden state helper.
        // ---------------------------------------------------------------------
        State candidate_state{
            committed[0],
            committed[1],
            committed[2],
            committed[3],
            committed[4],
            committed[5],
            committed[6]
        };

        Mat6 candidate_tangent;

        if constexpr (Finite) {
            // Finite-strain candidate: Green-Lagrange strain -> PK2 stress.
            const FiniteResponse response = integrate_finite_strain(
                volume_strain,
                candidate_state,
                shear,
                bulk,
                yield_curve
            );

            stress = response.stress;
            candidate_tangent = tangent_finite(
                response,
                shear,
                bulk,
                yield_curve
            );
        } else {
            // Small-strain candidate: infinitesimal strain -> Cauchy stress.
            const SmallResponse response = integrate_small_strain(
                volume_strain,
                candidate_state,
                shear,
                bulk,
                yield_curve
            );

            stress = response.stress;
            candidate_tangent = tangent_small(
                response,
                shear,
                bulk,
                yield_curve
            );
        }

        // ---------------------------------------------------------------------
        // Plane-stress residual and exact Newton derivative.
        //
        //     r     = S33
        //     dr/de33 = C3333.
        // ---------------------------------------------------------------------
        const Precision residual = stress[VolumeStress::Component::ZZ];

        if (std::abs(residual) <= tolerance) {
            converged = true;

            converged_state[0] = candidate_state[0];
            converged_state[1] = candidate_state[1];
            converged_state[2] = candidate_state[2];
            converged_state[3] = candidate_state[3];
            converged_state[4] = candidate_state[4];
            converged_state[5] = candidate_state[5];
            converged_state[6] = candidate_state[6];

            converged_tangent = candidate_tangent;
            break;
        }

        const Precision derivative = candidate_tangent(2, 2);
        logging::error(std::isfinite(derivative)
                       && std::abs(derivative)
                          > Precision(100) * std::numeric_limits<Precision>::epsilon(),
            "J2: singular shell plane-stress Newton derivative");

        const Precision delta = -residual / derivative;

        if constexpr (Finite) {
            // -----------------------------------------------------------------
            // Keep the Newton correction inside the admissible Green-Lagrange
            // domain. Only kinematics are checked here; the constitutive residual
            // itself is evaluated at the next Newton iteration.
            // -----------------------------------------------------------------
            Precision scale = Precision(1);
            bool accepted = false;

            while (scale >= Precision(1e-8)) {
                VolumeStrain candidate_strain = volume_strain;
                candidate_strain.voigt()(2) = e33 + scale * delta;

                const Mat3 candidate_C = Mat3::Identity()
                                       + Precision(2) * candidate_strain.tensor();
                Eigen::SelfAdjointEigenSolver<Mat3> solver(
                    Precision(0.5) * (candidate_C + candidate_C.transpose()),
                    Eigen::EigenvaluesOnly
                );

                if (solver.info() == Eigen::Success
                    && solver.eigenvalues().minCoeff() > Precision(1e-12)) {
                    e33 += scale * delta;
                    accepted = true;
                    break;
                }

                scale *= Precision(0.5);
            }

            logging::error(accepted,
                "J2: shell plane-stress Newton could not find admissible step");
        } else {
            e33 += delta;
        }
    }

    logging::error(converged,
        "J2: shell plane-stress constitutive solve did not converge");

    return stress;
}

/**
 * @brief Solves the two transverse strains required by a uniaxial stress state.
 *
 * Axial truss material response prescribes only e11. The transverse strains e22
 * and e33 are determined from the two stress conditions
 *
 *     S22 = 0
 *     S33 = 0.
 *
 * A two-dimensional Newton iteration is used. Its Jacobian is the transverse
 * block of the same analytically consistent three-dimensional J2 tangent,
 *
 *            [ C2222  C2233 ]
 *     J  =   [               ]
 *            [ C3322  C3333 ].
 *
 * After convergence the scalar axial tangent is obtained outside this routine by
 * Schur condensation of the same three-dimensional tangent.
 *
 * @tparam Finite false for linearized/Cauchy response, true for
 *                Green-Lagrange/PK2 response.
 * @param axial_strain Prescribed axial strain e11.
 * @param committed Seven-component committed J2 state.
 * @param converged_state Receives the state associated with the converged
 *                        transverse strains.
 * @param youngs Young's modulus E, used to scale the local tolerance.
 * @param poisson Poisson ratio nu, used for the elastic initial guess.
 * @param shear Elastic shear modulus G.
 * @param bulk Elastic bulk modulus K.
 * @param yield_curve Piecewise-linear isotropic hardening law.
 * @param converged_tangent Receives the full three-dimensional consistent tangent.
 * @return Converged full three-dimensional stress.
 */
template<bool Finite>
VolumeStress solve_axial_stress(Precision axial_strain,
                                const State& committed,
                                State& converged_state,
                                Precision youngs,
                                Precision poisson,
                                Precision shear,
                                Precision bulk,
                                const YieldCurve& yield_curve,
                                Mat6& converged_tangent) {
    // -------------------------------------------------------------------------
    // Elastic uniaxial-stress estimate for the two free transverse strains.
    // -------------------------------------------------------------------------
    VolumeStrain volume_strain;
    volume_strain.voigt().setZero();
    volume_strain.voigt()(0) = axial_strain;
    volume_strain.voigt()(1) = -poisson * axial_strain;
    volume_strain.voigt()(2) = -poisson * axial_strain;

    if constexpr (Finite) {
        volume_strain.voigt()(1) = std::max(
            volume_strain.voigt()(1), Precision(-0.49)
        );
        volume_strain.voigt()(2) = std::max(
            volume_strain.voigt()(2), Precision(-0.49)
        );
    }

    const Precision tolerance = stress_tolerance(
        youngs, initial_yield_stress(yield_curve)
    );

    VolumeStress stress;
    bool converged = false;

    for (Index iteration = 0; iteration < 30; ++iteration) {
        if constexpr (Finite) {
            // Check the finite-strain kinematic domain C = I + 2 E > 0.
            const Mat3 C = Mat3::Identity()
                         + Precision(2) * volume_strain.tensor();
            Eigen::SelfAdjointEigenSolver<Mat3> solver(
                Precision(0.5) * (C + C.transpose()),
                Eigen::EigenvaluesOnly
            );

            logging::error(solver.info() == Eigen::Success
                           && solver.eigenvalues().minCoeff() > Precision(1e-12),
                "J2: axial constitutive iterate left admissible Green-Lagrange domain");
        }

        // ---------------------------------------------------------------------
        // Start every transverse Newton candidate from the same committed
        // material history. Only the converged transverse solution may become the
        // caller's trial material state.
        // ---------------------------------------------------------------------
        State candidate_state{
            committed[0],
            committed[1],
            committed[2],
            committed[3],
            committed[4],
            committed[5],
            committed[6]
        };

        Mat6 candidate_tangent;

        if constexpr (Finite) {
            const FiniteResponse response = integrate_finite_strain(
                volume_strain,
                candidate_state,
                shear,
                bulk,
                yield_curve
            );

            stress = response.stress;
            candidate_tangent = tangent_finite(
                response,
                shear,
                bulk,
                yield_curve
            );
        } else {
            const SmallResponse response = integrate_small_strain(
                volume_strain,
                candidate_state,
                shear,
                bulk,
                yield_curve
            );

            stress = response.stress;
            candidate_tangent = tangent_small(
                response,
                shear,
                bulk,
                yield_curve
            );
        }

        // ---------------------------------------------------------------------
        // Residual of the uniaxial stress reduction.
        // ---------------------------------------------------------------------
        Eigen::Matrix<Precision, 2, 1> residual;
        residual << stress[VolumeStress::Component::YY],
                    stress[VolumeStress::Component::ZZ];

        if (residual.cwiseAbs().maxCoeff() <= tolerance) {
            converged = true;

            converged_state[0] = candidate_state[0];
            converged_state[1] = candidate_state[1];
            converged_state[2] = candidate_state[2];
            converged_state[3] = candidate_state[3];
            converged_state[4] = candidate_state[4];
            converged_state[5] = candidate_state[5];
            converged_state[6] = candidate_state[6];

            converged_tangent = candidate_tangent;
            break;
        }

        // ---------------------------------------------------------------------
        // Exact transverse Jacobian extracted from the consistent 3D tangent.
        // ---------------------------------------------------------------------
        Eigen::Matrix<Precision, 2, 2> jacobian;
        jacobian << candidate_tangent(1, 1), candidate_tangent(1, 2),
                    candidate_tangent(2, 1), candidate_tangent(2, 2);

        Eigen::FullPivLU<Eigen::Matrix<Precision, 2, 2>> lu(jacobian);
        logging::error(lu.isInvertible(),
            "J2: singular transverse Jacobian in axial stress reduction");

        const Eigen::Matrix<Precision, 2, 1> delta = lu.solve(-residual);

        if constexpr (Finite) {
            // Keep the transverse Newton correction in the admissible metric
            // domain before accepting the next constitutive candidate.
            Precision scale = Precision(1);
            bool accepted = false;

            while (scale >= Precision(1e-8)) {
                VolumeStrain candidate_strain = volume_strain;
                candidate_strain.voigt()(1) += scale * delta(0);
                candidate_strain.voigt()(2) += scale * delta(1);

                const Mat3 candidate_C = Mat3::Identity()
                                       + Precision(2) * candidate_strain.tensor();
                Eigen::SelfAdjointEigenSolver<Mat3> solver(
                    Precision(0.5) * (candidate_C + candidate_C.transpose()),
                    Eigen::EigenvaluesOnly
                );

                if (solver.info() == Eigen::Success
                    && solver.eigenvalues().minCoeff() > Precision(1e-12)) {
                    volume_strain = candidate_strain;
                    accepted = true;
                    break;
                }

                scale *= Precision(0.5);
            }

            logging::error(accepted,
                "J2: axial constitutive Newton could not find admissible step");
        } else {
            volume_strain.voigt()(1) += delta(0);
            volume_strain.voigt()(2) += delta(1);
        }
    }

    logging::error(converged,
        "J2: axial transverse-stress constitutive solve did not converge");

    return stress;
}
