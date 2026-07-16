#pragma once

#include "elasticity.h"

namespace fem::material {

struct ABDElasticity : Elasticity {
    Mat6 abd   = Mat6::Zero();
    Mat2 shear = Mat2::Zero();

    ABDElasticity(const Mat6& abd_in, const Mat2& shear_in);

    bool supports_shell_resultants() const override;

    void evaluate(const ShellGeneralizedStrain& strain,
                  Precision                     thickness,
                  const Precision*              state_old,
                  Precision*                    state_new,
                  ShellStressResultants&        resultants,
                  Mat8&                         tangent) const override;
};

} // namespace fem::material
