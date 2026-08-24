/**
 * @file reference.h
 * @brief Resolves reader-facing instance-qualified entity identifiers.
 *
 * Input decks address semantic entities before and after model compilation with
 * the same compact notation. An unqualified integer refers to the implicit
 * default instance, while `INSTANCE.ID` addresses the corresponding local entity
 * of a named instance. These helpers parse that scalar entity notation and map
 * it through the compiled Model identifier maps where required.
 *
 * Named set composition is intentionally not handled here. Adding node, element
 * or surface references to named regions belongs to `Model::add_nodes_to_set()`,
 * `Model::add_elements_to_set()` and `Model::add_surfaces_to_set()`, so parser
 * formats share the same model-level set semantics without a generic reader-side
 * container helper.
 *
 * @see model::Model
 *
 * @author Finn Eggers
 * @date 24.08.2026
 */

#pragma once

#include "../../core/logging.h"
#include "../../core/types_num.h"
#include "../../model/model.h"

#include <charconv>
#include <string>
#include <system_error>
#include <utility>

namespace fem::io::reader {

inline std::string qualify_reference(const std::string& target, const std::string& instance) {
    if (instance.empty() || target.find('.') != std::string::npos) {
        return target;
    }
    return instance + "." + target;
}

inline ID parse_local_id(const std::string& token, const std::string& context) {
    ID id{};
    const char* begin = token.data();
    const char* end   = begin + token.size();
    const auto [ptr, ec] = std::from_chars(begin, end, id);

    logging::error(ec == std::errc{} && ptr == end,
        context, ": '", token, "' is not a valid integer identifier");
    return id;
}

inline std::pair<std::string, ID> split_entity_reference(const std::string& target) {
    const auto separator = target.rfind('.');
    if (separator == std::string::npos) {
        return {model::Model::DEFAULT_INSTANCE_NAME, parse_local_id(target, "Entity reference")};
    }

    const std::string instance = target.substr(0, separator);
    const std::string local    = target.substr(separator + 1);

    logging::error(!instance.empty() && !local.empty(),
        "Entity reference '", target, "' must use INSTANCE.ID");
    return {instance, parse_local_id(local, "Entity reference")};
}

inline ID compiled_node_id(const model::Model& model, const std::string& target) {
    const auto [instance, id] = split_entity_reference(target);
    return model.compiled_node_id(id, instance);
}

inline ID compiled_element_id(const model::Model& model, const std::string& target) {
    const auto [instance, id] = split_entity_reference(target);
    return model.compiled_element_id(id, instance);
}

inline ID compiled_surface_id(const model::Model& model, const std::string& target) {
    const auto [instance, id] = split_entity_reference(target);
    return model.compiled_surface_id(id, instance);
}

inline ID compiled_line_id(const model::Model& model, const std::string& target) {
    const auto [instance, id] = split_entity_reference(target);
    return model.compiled_line_id(id, instance);
}

} // namespace fem::io::reader
