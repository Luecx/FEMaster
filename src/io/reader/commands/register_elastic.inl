/**
 * @file register_elastic.inl
 * @brief Registers supported linear-elastic material definitions.
 *
 * The command maps isotropic elasticity, generalized isotropic elasticity,
 * orthotropic engineering constants and orthotropic stiffness coefficients to
 * FEMaster elasticity models. `TYPE=ORTHOTROPIC` converts the supplied symmetric
 * normal stiffness block to engineering constants before constructing
 * `OrthotropicElasticity`.
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

#include <Eigen/LU>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../material/generalised_isotropic_elasticity.h"
#include "../../../material/isotropic_elasticity.h"
#include "../../../material/orthotropic_elasticity.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

/**
 * Registers linear-elastic properties for the active material.
 *
 * Supported forms are isotropic `E, nu`, generalized isotropic `E, nu, G`,
 * orthotropic engineering constants
 * `E1, E2, E3, nu12, nu13, nu23, G12, G13, G23`, and Abaqus orthotropic
 * stiffness coefficients
 * `D1111, D1122, D2222, D1133, D2233, D3333, D1212, D1313, D2323`.
 *
 * @param registry DSL registry receiving the elastic command.
 * @param model FEMaster model containing the active material.
 */
inline void register_elastic(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ELASTIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign linear-elastic properties to the active material.");

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
                    logging::error(material != nullptr,
                        "ELASTIC requires an active material context");
                    material->set_elasticity<fem::material::IsotropicElasticity>(E, nu);
                })
            )
        );

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
                    logging::error(material != nullptr,
                        "ELASTIC requires an active material context");
                    material->set_elasticity<fem::material::GeneralisedIsotropicElasticity>(E, nu, G);
                })
            )
        );

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
                    logging::error(material != nullptr,
                        "ELASTIC requires an active material context");

                    material->set_elasticity<fem::material::OrthotropicElasticity>(
                        data[0], data[1], data[2],
                        data[3], data[4], data[5],
                        data[6], data[7], data[8]
                    );
                })
            )
        );

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
                    logging::error(material != nullptr,
                        "ELASTIC requires an active material context");

                    fem::Mat3 normal_stiffness;
                    normal_stiffness << data[0], data[1], data[3],
                                        data[1], data[2], data[4],
                                        data[3], data[4], data[5];

                    const fem::Precision scale = normal_stiffness.cwiseAbs().maxCoeff();
                    logging::error(scale > fem::Precision(0),
                        "ELASTIC TYPE=ORTHOTROPIC has a singular normal stiffness block");

                    const fem::Mat3 normalized = normal_stiffness / scale;
                    logging::error(std::abs(normalized.determinant()) >
                                std::numeric_limits<fem::Precision>::epsilon(),
                        "ELASTIC TYPE=ORTHOTROPIC has a singular normal stiffness block");

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
