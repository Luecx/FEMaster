    // Finite-strain response with the identical constant material operator,
    // interpreted as the mapping from Green-Lagrange strain to PK2 stress.
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6*                            tangent = nullptr) const override;

    // Linearized shell plane-stress response. The in-plane normal block uses
    // E and nu; in-plane and transverse engineering shear terms use G.
    // Finite-strain shell response returning PK2 components. The reduced
    // material derivative is written only when requested.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5*                                   tangent = nullptr) const override;

private:
    // Build the generalized-isotropic in-plane plane-stress operator.
    [[nodiscard]] Mat3 plane_stress_tangent() const;

    // Embed in-plane and transverse shear terms in shell material ordering.
    [[nodiscard]] Mat5 shell_material_tangent() const;

    // Build the three-dimensional normal-coupling and independent-shear tangent.
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
