/**
 * @file register_tload.inl
 * @brief Registers thermal loads driven by nodal temperature fields.
 *
 * `TLOAD` associates a compiled nodal temperature field and reference
 * temperature with a named load collector. The resulting `bc::TLoad` provides
 * the temperature difference consumed by elements whose materials define a
 * thermal-expansion coefficient.
 *
 * This registration does not compute thermal strains or equivalent forces;
 * those operations remain coupled to element geometry and constitutive data
 * during analysis.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../../bc/load_t.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_tload(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("TLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").required()
                .key("TEMPERATUREFIELD").required()
                .key("REFERENCETEMPERATURE").required()
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->load_cols.activate(keys.raw("LOAD_COLLECTOR"));

            const std::string field_name = keys.raw("TEMPERATUREFIELD");
            const auto field = model._data->get_field(field_name);
            logging::error(field != nullptr,
                "TLOAD: temperature field ", field_name, " does not exist");
            logging::error(field->domain == model::FieldDomain::NODE,
                "TLOAD: temperature field ", field_name, " must use NODE domain");
            logging::error(field->components == 1,
                "TLOAD: temperature field ", field_name, " must have one component");

            auto load = std::make_shared<bc::TLoad>();
            load->temp_field_ = field;
            load->ref_temp_   = keys.get<fem::Precision>("REFERENCETEMPERATURE");
            model.add_load(std::move(load));
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
