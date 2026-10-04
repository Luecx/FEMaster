    // Finite-strain shell response returning PK2 components work-conjugate to
    // the five supplied Green-Lagrange strain components. State remains unchanged.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5*                                   tangent = nullptr) const override;

private:
    // Build the in-plane orthotropic plane-stress tangent ordered [11,22,12].
    [[nodiscard]] Mat3 plane_stress_tangent() const;

    // Embed the in-plane block and directional transverse shear moduli into the
    // five-component shell material ordering.
    [[nodiscard]] Mat5 shell_material_tangent() const;

    // Invert the symmetric engineering compliance to obtain the full
    // three-dimensional tangent in local material coordinates.
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
