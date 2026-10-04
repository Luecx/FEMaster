bool OrthotropicElasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

/**
 * Builds the orthotropic in-plane plane-stress tangent.
 *
 * The reciprocal ratio `nu21 = nu12 E2/E1` enforces symmetry of the normal
 * block. The engineering 1-2 shear term remains uncoupled and uses `G12`.
 *
 * @return Constant tangent ordered `[11,22,12]`.
 */
Mat3 OrthotropicElasticity::plane_stress_tangent() const {
    const Precision nu21  = nu12 * E2 / E1;
    const Precision denom = Precision(1) - nu12 * nu21;

    Mat3 tangent;
    tangent << E1 / denom,         nu12 * E2 / denom, Precision(0),
               nu21 * E1 / denom, E2 / denom,        Precision(0),
               Precision(0),       Precision(0),       G12;
    return tangent;
}

/**
 * Embeds the plane-stress block and directional transverse shear moduli into
 * the five-component shell material tangent.
 *
 * @return Tangent ordered `[11,22,12,13,23]`.
 */
Mat5 OrthotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = G13;
    tangent(4, 4) = G23;
    return tangent;
}

/**
 * Constructs the full three-dimensional orthotropic tangent.
 *
 * The engineering compliance is assembled directly from the three major
 * Poisson ratios. Its off-diagonal terms are written in their symmetric form,
 * so the reciprocal minor ratios never become independent stored parameters.
 * FEMaster's volume Voigt ordering places shear components as `[23,13,12]`.
 *
 * @return Constant six-by-six engineering-Voigt tangent.
 */
Mat6 OrthotropicElasticity::volume_tangent() const {
    Mat6 compliance;
    compliance <<
        Precision(1) / E1, -nu12 / E1,         -nu13 / E1,         Precision(0),       Precision(0),       Precision(0),
        -nu12 / E1,         Precision(1) / E2, -nu23 / E2,         Precision(0),       Precision(0),       Precision(0),
        -nu13 / E1,         -nu23 / E2,         Precision(1) / E3, Precision(0),       Precision(0),       Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(1) / G23, Precision(0),       Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(1) / G13, Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(1) / G12;
    return compliance.inverse();
}

/**
 * Evaluates orthotropic PK2 stress from Green-Lagrange strain.
 *
 * @param strain Green-Lagrange engineering strain vector.
 * @param old_state Unused input material-point state row.
 * @param new_state Unused output material-point state row.
 * @param stress Second Piola-Kirchhoff stress in material coordinates.
 * @param tangent Optional material derivative `dS/dE`.
 */
void OrthotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     const Precision*                 old_state,
                                     Precision*                       new_state,
                                     VolumeStressPK2&                 stress,
                                     Mat6*                            tangent) const {
    (void) old_state;
    (void) new_state;

    const Mat6 material_tangent = volume_tangent();
    stress.voigt() = material_tangent * strain.voigt();

    if (tangent != nullptr) {
        *tangent = material_tangent;
    }
}

/**
 * Evaluates orthotropic shell PK2 stress from Green-Lagrange strain.
 *
 * @param strain Five-component Green-Lagrange material strain.
 * @param old_state Unused input material-point state row.
 * @param new_state Unused output material-point state row.
 * @param stress Shell second Piola-Kirchhoff stress.
 * @param tangent Optional reduced material derivative.
 */
void OrthotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     const Precision*                        old_state,
                                     Precision*                              new_state,
                                     ShellMaterialStressPK2&                 stress,
                                     Mat5*                                   tangent) const {
    (void) old_state;
    (void) new_state;

    const Mat5 material_tangent = shell_material_tangent();
    stress.values() = material_tangent * strain.values();

    if (tangent != nullptr) {
        *tangent = material_tangent;
    }
}

} // namespace fem::material
