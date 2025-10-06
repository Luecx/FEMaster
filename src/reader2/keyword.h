/**
 * @file keyword.h
 * @brief Strongly-typed representation of a star-keyword line (e.g. "*MATERIAL, NAME=STEEL").
 *
 * A Keyword carries the command name, the scope at which it was read, and a validated
 * key->value map. Construction centralizes parsing from a Line and KeyRules validation.
 */
#pragma once

#include <string>
#include <unordered_map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include "types.h"
#include "line.h"
#include "key_rules.h"

namespace fem::reader2 {

/// Represents a parsed keyword line within the current parser scope.
class Keyword {
public:
    using Map = std::unordered_map<std::string, std::string>;

    /**
     * @brief Construct a Keyword from a parsed KEYWORD line, current scope, and optional rules.
     * @throws std::runtime_error if @p L is not a KEYWORD line or key validation fails.
     */
    static Keyword from_line(const Line& L,
                             const Scope& current_scope,
                             const std::optional<KeyRules>& rules);

    /// Current scope (e.g. "ROOT", "MATERIAL", ...).
    const Scope& scope()  const { return scope_; }
    /// Uppercased command name (e.g. "MATERIAL", "ELASTIC", ...).
    const std::string& name()   const { return name_;  }
    /// Validated key map (uppercased keys).
    const Map& kv()      const { return kv_; }

    /// True if a key exists.
    bool has(const std::string& k) const { return kv_.find(k) != kv_.end(); }

    /// Typed getter with default.
    template<class T>
    T get(const std::string& k, const T& def) const {
        auto it = kv_.find(k);
        if (it == kv_.end()) return def;
        return convert<T>(it->second);
    }

    /// Typed getter, throws if missing.
    template<class T>
    T require(const std::string& k) const {
        auto it = kv_.find(k);
        if (it == kv_.end()) throw std::runtime_error("Missing required key: " + k);
        return convert<T>(it->second);
    }

private:
    Scope        scope_;
    std::string  name_;
    Map          kv_;

    Keyword(Scope s, std::string n, Map m)
        : scope_(std::move(s)), name_(std::move(n)), kv_(std::move(m)) {}

    template<class T>
    static T convert(const std::string& s) {
        if constexpr (std::is_same_v<T, std::string>) {
            return s;
        } else {
            std::istringstream iss(s);
            T v{};
            iss >> v;
            if (iss.fail()) throw std::runtime_error("Key conversion failed for value '" + s + "'");
            return v;
        }
    }
};

} // namespace fem::reader2
