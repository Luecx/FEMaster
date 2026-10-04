    // Total-Lagrangian three-dimensional response in the reference material
    // basis. Green-Lagrange strain, PK2 stress and dS/dE remain work-conjugate;
    // the owning section handles transformations to and from global coordinates.
    virtual void evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                 old_state,
                          Precision*                       new_state,
                          VolumeStressPK2&                 stress,
                          Mat6*                            tangent = nullptr) const;

    // Generalized beam response. The section-defined six-component strain and
    // resultant ordering is preserved. The optional tangent is their consistent
    // local derivative.
    virtual void evaluate(const BeamGeneralizedStrain& strain,
                          const Precision*             old_state,
                          Precision*                   new_state,
                          BeamStressResultants&        resultants,
                          Mat6*                        tangent = nullptr) const;

    // Infinitesimal shell material response at one physical thickness point.
    // The five strain components exclude thickness-normal strain; stress is
    // Cauchy stress under the material's plane-stress reduction.
    // Finite-strain shell material response at one physical thickness point.
    // The five-component Green-Lagrange input returns work-conjugate PK2 stress
    // and, when requested, the consistently reduced tangent under S33 = 0.
    virtual void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                        old_state,
                          Precision*                              new_state,
                          ShellMaterialStressPK2&                 stress,
                          Mat5*                                   tangent = nullptr) const;

    // Runtime access to concrete constitutive implementations
    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }

    template<typename T>
    const T* as() const {
        return dynamic_cast<const T*>(this);
    }
};

using ElasticityPtr = Elasticity::Ptr;

} // namespace material
} // namespace fem
