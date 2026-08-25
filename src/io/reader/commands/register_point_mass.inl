/**
 * @file register_point_mass.inl
 * @brief Registers concentrated FEMaster point masses on compiled node sets.
 *
 * The public `POINTMASS, NSET=...` syntax is unchanged. The command executes on
 * the compiled model, resolves the named node region and creates one auxiliary
 * `PointElement` per target node. All physical coefficients are stored in one
 * shared `PointMassSection`; no PointMass feature object is created.
 *
 * The generated point elements deliberately stay outside the frozen dense
 * element namespace because the command runs after `Model::compile()`. They
 * participate in DOF, stiffness, mass and inertia assembly without changing
 * ELSETs or element-domain field enumeration.
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../model/element/point.h"
#include "../../../model/model.h"
#include "../../../section/section_point_mass.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers the existing NSET-based FEMaster point-mass command.
 *
 * The data layout remains unchanged: scalar translational mass, three rotary
 * inertias, three translational ground stiffnesses and three rotational ground
 * stiffnesses.
 *
 * @param registry Parser registry receiving the command.
 * @param model Compiled model receiving auxiliary point elements.
 */
inline void register_point_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("POINTMASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign concentrated mass, inertia and ground stiffness to a compiled node set.");

        auto nset = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET").required().doc("Target compiled node set")
        );

        command.on_enter([nset](const fem::io::dsl::Keys& keys) {
            *nset = keys.raw("NSET");
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is("ROOT"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::Precision>().name("MASS").desc("Point mass")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("INERTIA").desc("Rotary inertia Ix,Iy,Iz")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("SPRING").desc("Translational spring constants")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("ROTSPRING").desc("Rotational spring constants")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, nset](
                          fem::Precision mass,
                          const std::array<fem::Precision, 3>& inertia_data,
                          const std::array<fem::Precision, 3>& spring_data,
                          const std::array<fem::Precision, 3>& rotary_data) {
                    logging::error(model._data != nullptr && model._data->compiled,
                        "POINTMASS: command requires a compiled model");
                    logging::error(model._data->node_sets.has(*nset),
                        "POINTMASS: node set ", *nset, " does not exist");

                    const auto nodes = model._data->node_sets.get(*nset);
                    logging::error(nodes != nullptr,
                        "POINTMASS: node set ", *nset, " is not initialized");

                    const fem::Vec3 inertia{inertia_data[0], inertia_data[1], inertia_data[2]};
                    const fem::Vec3 spring {spring_data[0],  spring_data[1],  spring_data[2]};
                    const fem::Vec3 rotary {rotary_data[0],  rotary_data[1],  rotary_data[2]};

                    auto section_region = std::make_shared<model::ElementRegion>(
                        "__POINTMASS_" + std::to_string(model._data->point_elements.size()));
                    auto section = std::make_shared<PointMassSection>(
                        section_region, mass, inertia, spring, rotary);

                    for (const ID node : *nodes) {
                        const ID id = -1 - static_cast<ID>(model._data->point_elements.size());
                        section_region->add(id);

                        auto point = std::make_shared<model::PointElement>(
                            id, std::array<ID, model::PointElement::N>{node});
                        point->_model_data = model._data.get();
                        point->set_section(section);
                        model._data->point_elements.push_back(std::move(point));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
