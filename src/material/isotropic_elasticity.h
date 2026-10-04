    // Finite-strain response using the same constant Hooke operator. Input is
    // Green-Lagrange strain and output is second Piola-Kirchhoff stress.
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6*                            tangent = nullptr) const override;

    // Linearized five-component shell response with an in-plane plane-stress
    // block and transverse shear moduli.
    // Green-Lagrange five-component shell response with PK2 output. The material
    // state remains unchanged and the reduced tangent is optional.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5*                                   tangent = nullptr) const override;

private:
    // Build the in-plane plane-stress operator ordered as [11,22,12].
    [[nodiscard]] Mat3 plane_stress_tangent() const;

    // Embed the plane-stress block and transverse shear moduli into the
    // five-component shell material ordering [11,22,12,13,23].
    [[nodiscard]] Mat5 shell_material_tangent() const;

    // Build the full isotropic three-dimensional engineering-Voigt tangent.
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
