/**
 * @file reference.h
 * @brief Provides reader-side normalization of entity-reference tokens.
 *
 * Input commands may provide an entity token together with a separate
 * `INSTANCE` keyword. `qualify_reference()` combines those syntactic components
 * without replacing an Instance already present in the token. Commands that
 * construct Part-local topology can use `parse_local_id()` for integer tokens
 * that are not handled by a model-level region operation.
 *
 * Interpretation of `ID` and `INSTANCE.ID`, mapping through compiled Instance
 * data and named-region composition belong to `model::Model`. This file retains
 * only transformations that are specific to reader grammar.
 *
 * @author Finn Eggers
 * @date 24.08.2026
 */

#pragma once

#include "../../core/logging.h"
#include "../../core/types_num.h"

#include <charconv>
#include <string>
#include <system_error>

namespace fem::io::reader {

// Apply a separate Instance keyword only to an otherwise unqualified token
inline std::string qualify_reference(const std::string& target, const std::string& instance) {
    if (instance.empty() || target.find('.') != std::string::npos) {
        return target;
    }
    return instance + "." + target;
}

// Parse a complete Part-local integer token for commands that construct topology
inline ID parse_local_id(const std::string& token, const std::string& context) {
    ID id{};
    const char* begin = token.data();
    const char* end   = begin + token.size();
    const auto [ptr, ec] = std::from_chars(begin, end, id);

    logging::error(ec == std::errc{} && ptr == end,
        context, ": '", token, "' is not a valid integer identifier");
    return id;
}

} // namespace fem::io::reader
