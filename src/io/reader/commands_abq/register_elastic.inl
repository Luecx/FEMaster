/**
 * @file register_elastic.inl
 * @brief Registers Abaqus *ELASTIC material definitions supported by FEMaster.
 *
 * Abaqus and the native FEMaster reader use different meanings for
 * `TYPE=ORTHOTROPIC`, so Abaqus elasticity is parsed independently. Isotropic
 * elasticity, engineering-constant orthotropy and Abaqus orthotropic stiffness
 * coefficients are converted into the existing FEMaster material models.
 *
 * Abaqus lamina elasticity, transverse isotropy and fully anisotropic elasticity
 * remain unsupported until an exact FEMaster representation is added
 * deliberately.
 *
 * @see material::IsotropicElasticity
 * @see material::OrthotropicElasticity
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/LU>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../material/isotropic_elasticity.h"
#include "../../../material/orthotropic_elasticity.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus linear-elastic material data for the active material.
 *
 * `TYPE=ISOTROPIC` is the Abaqus default and reads Young's modulus and Poisson's
 * ratio. `TYPE=ENGINEERING CONSTANTS` reads the Abaqus ordering
 *
 * `E1, E2, E3, nu12, nu13, nu23, G12, G13, G23`.
 *
 * `TYPE=ORTHOTROPIC` reads the nine independent Abaqus stiffness coefficients.
 * The symmetric normal 3x3 block is inverted to recover engineering compliance,
 * from which the directional Young's moduli and Poisson ratios are obtained.
 * The three uncoupled shear coefficients map directly to the corresponding
 * engineering shear moduli.
 *
 * Temperature and field-variable dependencies are intentionally unsupported in
 * this initial reader and therefore cannot be supplied as additional data.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model containing the active material.
 */
inline void register_elastic(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ELASTIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign supported Abaqus linear-elastic properties to the active material.");

        // Accept only Abaqus elasticity forms that map exactly onto an existing
        // FEMaster constitutive implementation
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE")
                    .optional("ISOTROPIC")
                    .allowed({"ISOTROPIC", "ENGINEERINGCONSTANTS", "ORTHOTROPIC"})
                    .doc("Abaqus elastic material form")
        );

        // Map the Abaqus isotropic pair directly to FEMaster isotropic Hooke
        // elasticity
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ISOTROPIC"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("E").desc("Young's modulus")
                    .one<fem::Precision>().name("NU").desc("Poisson ratio")
                )
                .bind([&model](fem::Precision E, fem::Precision nu) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("ELASTIC requires an active material context");
                    }
                    material->set_elasticity<fem::material::IsotropicElasticity>(E, nu);
                })
            )
        );

        // Reorder Abaqus engineering constants into the FEMaster orthotropic
        // constructor convention and derive the reciprocal Poisson ratio nu31
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ENGINEERINGCONSTANTS"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(2))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>().name("DATA")
                        .desc("E1,E2,E3,nu12,nu13,nu23,G12,G13,G23")
                )
                .bind([&model](const std::array<fem::Precision, 9>& data) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("ELASTIC requires an active material context");
                    }

                    const fem::Precision E1   = data[0];
                    const fem::Precision E2   = data[1];
                    const fem::Precision E3   = data[2];
                    const fem::Precision nu12 = data[3];
                    const fem::Precision nu13 = data[4];
                    const fem::Precision nu23 = data[5];
                    const fem::Precision G12  = data[6];
                    const fem::Precision G13  = data[7];
                    const fem::Precision G23  = data[8];
                    const fem::Precision nu31 = nu13 * E3 / E1;

                    material->set_elasticity<fem::material::OrthotropicElasticity>(
                        E1, E2, E3, G23, G13, G12, nu23, nu31, nu12
                    );
                })
            )
        );

        // Recover engineering constants from the Abaqus orthotropic stiffness
        // coefficients before constructing the FEMaster orthotropic material
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ORTHOTROPIC"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(2))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>().name("DATA")
                        .desc("D1111,D1122,D2222,D1133,D2233,D3333,D1212,D1313,D2323")
                )
                .bind([&model](const std::array<fem::Precision, 9>& data) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("ELASTIC requires an active material context");
                    }

                    // Assemble and invert the normal stiffness block to recover
                    // the symmetric engineering compliance coefficients
                    fem::Mat3 normal_stiffness;
                    normal_stiffness << data[0], data[1], data[3],
                                        data[1], data[2], data[4],
                                        data[3], data[4], data[5];

                    const fem::Precision determinant = normal_stiffness.determinant();
                    if (std::abs(determinant) <= std::numeric_limits<fem::Precision>::epsilon()) {
                        throw std::runtime_error("ELASTIC TYPE=ORTHOTROPIC has a singular normal stiffness block");
                    }
                    const fem::Mat3 compliance = normal_stiffness.inverse();

                    const fem::Precision E1   = fem::Precision(1) / compliance(0, 0);
                    const fem::Precision E2   = fem::Precision(1) / compliance(1, 1);
                    const fem::Precision E3   = fem::Precision(1) / compliance(2, 2);
                    const fem::Precision nu12 = -compliance(0, 1) * E1;
                    const fem::Precision nu31 = -compliance(0, 2) * E3;
                    const fem::Precision nu23 = -compliance(1, 2) * E2;
                    const fem::Precision G12  = data[6];
                    const fem::Precision G13  = data[7];
                    const fem::Precision G23  = data[8];

                    material->set_elasticity<fem::material::OrthotropicElasticity>(
                        E1, E2, E3, G23, G13, G12, nu23, nu31, nu12
                    );
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
