// parser/keyword.cpp
/**
 * @file keyword.cpp
 * @brief Implementation of Keyword construction & validation.
 */
#include "keyword.h"

namespace fem::reader2 {

    Keyword Keyword::from_line(const Line& L,
                               const Scope& current_scope,
                               const std::optional<KeyRules>& rules)
    {
        if (L.type() != LineType::KEYWORD)
            throw std::runtime_error("Keyword::from_line: Line is not a KEYWORD");

        std::string name = L.command();
        Keyword::Map keys = L.keys(); // uppercased keys/values normalized by Line

        if (rules) {
            keys = rules->apply(keys);
        }

        return Keyword(current_scope, std::move(name), std::move(keys));
    }

} // namespace fem::reader2
