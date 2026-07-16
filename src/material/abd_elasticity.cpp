#include "abd_elasticity.h"

#include "strain/shell_generalized_strain.h"
#include "stress/shell_stress_resultants.h"

namespace fem::material {

ABDElasticity::ABDElasticity(const Mat6& abd_in, const Mat2& shear_in)
    : abd  (abd_in),
      shear(shear_in) {}

bool ABDElasticity::supports_shell_resultants() const {
    return true;
}

void ABDElasticity::evaluate(const ShellGeneralizedStrain& strain,
                             Precision                     thickness,
                             const Precision*              state_old,
                             Precision*                    state_new,
                             ShellStressResultants&        resultants,
                             Mat8&                         tangent) const {
    (void) thickness;
    (void) state_old;
    (void) state_new;

    tangent.setZero();
    tangent.template block<6, 6>(0, 0) = abd;
    tangent.template block<2, 2>(6, 6) = shear;
    resultants.values() = tangent * strain.values();
}

} // namespace fem::material
