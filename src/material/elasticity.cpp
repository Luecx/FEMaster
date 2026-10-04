bool Elasticity::supports_axial_green_lagrange() const {
    return false;
}

bool Elasticity::supports_volume_green_lagrange() const {
    return false;
}

bool Elasticity::supports_beam_resultants() const {
    return false;
}

bool Elasticity::supports_shell_integration_green_lagrange() const {
    return false;
}

Index Elasticity::state_size() const {
    return 0;
}

void Elasticity::initialize_state(Precision* state) const {
    (void) state;
}

void Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                          const Precision*                old_state,
                          Precision*                      new_state,
                          AxialStressPK2&                 stress,
                          Precision*                      tangent) const {
    (void) strain;
    (void) old_state;
    (void) new_state;
    (void) stress;
    (void) tangent;

    logging::error(false,
        "Elasticity model does not support Green-Lagrange axial evaluation");
}

void Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                 old_state,
                          Precision*                       new_state,
                          VolumeStressPK2&                 stress,
                          Mat6*                            tangent) const {
    (void) strain;
    (void) old_state;
    (void) new_state;
    (void) stress;
    (void) tangent;

    logging::error(false,
        "Elasticity model does not support Green-Lagrange volume evaluation");
}

void Elasticity::evaluate(const BeamGeneralizedStrain& strain,
                          const Precision*             old_state,
                          Precision*                   new_state,
                          BeamStressResultants&        resultants,
                          Mat6*                        tangent) const {
    (void) strain;
    (void) old_state;
    (void) new_state;
    (void) resultants;
    (void) tangent;

    logging::error(false,
        "Elasticity model does not support beam-resultant evaluation");
}

void Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                        old_state,
                          Precision*                              new_state,
                          ShellMaterialStressPK2&                 stress,
                          Mat5*                                   tangent) const {
    (void) strain;
    (void) old_state;
    (void) new_state;
    (void) stress;
    (void) tangent;

    logging::error(false,
        "Elasticity model does not support Green-Lagrange shell integration");
}

} // namespace fem::material
