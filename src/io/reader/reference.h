/**
 * @file reference.h
 * @brief Resolves reader-facing instance-qualified entity references.
 *
 * Input decks address semantic entities before and after model compilation with
 * the same compact notation. An unqualified integer refers to the implicit
 * default instance, while `INSTANCE.ID` addresses the corresponding local entity
 * of a named instance. Named set references may use the same optional instance
 * prefix and are resolved from the compiled `ModelData` registries.
 *
 * These helpers deliberately remain free functions. They do not introduce a
 * parallel reference-object hierarchy and keep compiled identifiers confined to
 * parser/model internals.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../core/logging.h"
#include "../../core/types_num.h"
#include "../../data/sets.h"
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

template<typename Region, typename Resolver>
inline void add_compiled_reference(model::Sets<Region>& sets,
                                   typename Region::Ptr destination,
                                   const std::string& target,
                                   const std::string& instance,
                                   Resolver&& resolver) {
    logging::error(destination != nullptr,
        "Assembly set destination is not initialized");

    const std::string qualified = qualify_reference(target, instance);
    if (sets.has(qualified)) {
        auto source = sets.get(qualified);
        logging::error(source != nullptr,
            "Assembly set reference ", qualified, " is not initialized");
        if (source != destination) {
            destination->add(*source);
        }
        return;
    }

    destination->add(resolver(qualified));
}

} // namespace fem::io::reader
