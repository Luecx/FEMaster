/**
 * @file register_elastic.inl
 * @brief Registers shared FEMaster and Abaqus *ELASTIC material definitions.
 *
 * The common elastic command follows Abaqus linear-elastic type semantics so
 * native and Abaqus input decks use one material parser. Isotropic elasticity,
 * orthotropic engineering constants and direct orthotropic stiffness
 * coefficients are mapped onto the existing FEMaster elasticity models.
 *
 * `TYPE=ENGINEERING CONSTANTS` is the canonical external representation for
 * `OrthotropicElasticity`. `TYPE=ORTHOTROPIC` accepts the Abaqus stiffness
 * coefficients and converts their normal stiffness block to engineering
 * constants before constructing the same material model.
 *
 * The existing generalized-isotropic FEMaster type remains available as a
 * native extension because it represents a distinct constitutive model already
 * used by FEMaster.
 *
 * @see material::IsotropicElasticity
 * @see material::GeneralisedIsotropicElasticity
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
#include "../../../material/generalised_isotropic_elasticity.h"
#include "../../../material/isotropic_elasticity.h"
#include "../../../material/orthotropic_elasticity.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

/**
 * Registers linear-elastic material properties for the active material.
 *
 * The supported standard forms are:
 *
 * - `TYPE=ISOTROPIC`: `E, nu`,
 * - `TYPE=ENGINEERING CONSTANTS`:
 *   `E1, E2, E3, nu12, nu13, nu23, G12, G13, G23`,
 * - `TYPE=ORTHOTROPIC`:
 *   `D1111, D1122, D2222, D1133, D2233, D3333, D1212, D1313, D2323`.
 *
 * Direct orthotropic stiffness coefficients are reduced to the canonical
 * engineering-constants representation by inverting the symmetric 3x3 normal
 * stiffness block. The three uncoupled shear coefficients map directly to
 * `G12`, `G13` and `G23`.
 *
 * The legacy FEMaster `GENERALISED ISOTROPIC` aliases remain supported as an
 * extension and map to `GeneralisedIsotropicElasticity` without affecting the
 * standard Abaqus type meanings.
 *
 * Temperature- and field-dependent elastic data are not supported by this
 * registration and therefore must not be appended to the material rows.
 *
 * @param registry DSL registry receiving the elastic command.
 * @param model FEMaster model containing the active material.
 */
inline void register_elastic(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ELASTIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign linear-elastic properties to the active material.");

        // Use Abaqus type meanings for the standard elastic forms while keeping
        // the existing generalized-isotropic FEMaster extension available
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE")
                    .optional("ISOTROPIC")
                    .allowed({
                        "ISO", "ISOTROPIC",
                        "GENERALISEDISOTROPIC", "GENERALISED_ISOTROPIC", "GENISO",
                        "ENGINEERINGCONSTANTS",
                        "ORTHO", "ORTHOTROPIC"
                    })
                    .doc("Elasticity formulation")
        );

        // Map the isotropic Young's modulus and Poisson ratio directly
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ISO", "ISOTROPIC"}))
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

        // Preserve the existing FEMaster generalized-isotropic extension
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals(
                "TYPE", {"GENERALISEDISOTROPIC", "GENERALISED_ISOTROPIC", "GENISO"}
            ))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("E").desc("Young's modulus")
                    .one<fem::Precision>().name("NU").desc("Poisson ratio")
                    .one<fem::Precision>().name("G").desc("Independent shear modulus")
                )
                .bind([&model](fem::Precision E, fem::Precision nu, fem::Precision G) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("ELASTIC requires an active material context");
                    }
                    material->set_elasticity<fem::material::GeneralisedIsotropicElasticity>(E, nu, G);
                })
            )
        );

        // Construct orthotropic elasticity directly from conventional
        // engineering constants in Abaqus ordering
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

                    material->set_elasticity<fem::material::OrthotropicElasticity>(
                        data[0], data[1], data[2],
                        data[3], data[4], data[5],
                        data[6], data[7], data[8]
                    );
                })
            )
        );

        // Convert Abaqus orthotropic stiffness coefficients to engineering
        // constants before constructing the canonical orthotropic material
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ORTHO", "ORTHOTROPIC"}))
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

                    // Assemble the symmetric normal stiffness block and test its
                    // normalized determinant before inversion
                    fem::Mat3 normal_stiffness;
                    normal_stiffness << data[0], data[1], data[3],
                                        data[1], data[2], data[4],
                                        data[3], data[4], data[5];

                    const fem::Precision scale = normal_stiffness.cwiseAbs().maxCoeff();
                    if (scale <= fem::Precision(0)) {
                        throw std::runtime_error("ELASTIC TYPE=ORTHOTROPIC has a singular normal stiffness block");
                    }

                    const fem::Mat3 normalized = normal_stiffness / scale;
                    if (std::abs(normalized.determinant()) <= std::numeric_limits<fem::Precision>::epsilon()) {
                        throw std::runtime_error("ELASTIC TYPE=ORTHOTROPIC has a singular normal stiffness block");
                    }

                    // Invert the normal block and recover the three major
                    // Poisson ratios in the canonical engineering convention
                    const fem::Mat3 compliance = normal_stiffness.inverse();

                    const fem::Precision E1   = fem::Precision(1) / compliance(0, 0);
                    const fem::Precision E2   = fem::Precision(1) / compliance(1, 1);
                    const fem::Precision E3   = fem::Precision(1) / compliance(2, 2);
                    const fem::Precision nu12 = -compliance(0, 1) * E1;
                    const fem::Precision nu13 = -compliance(0, 2) * E1;
                    const fem::Precision nu23 = -compliance(1, 2) * E2;
                    const fem::Precision G12  = data[6];
                    const fem::Precision G13  = data[7];
                    const fem::Precision G23  = data[8];

                    material->set_elasticity<fem::material::OrthotropicElasticity>(
                        E1, E2, E3, nu12, nu13, nu23, G12, G13, G23
                    );
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
