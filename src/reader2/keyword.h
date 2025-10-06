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
    const Scope& scope() const
    {
        return _scope;
    }
    /// Uppercased command name (e.g. "MATERIAL", "ELASTIC", ...).
    const std::string& name() const
    {
        return _name;
    }
    /// Validated key map (uppercased keys).
    const Map& kv() const
    {
        return _kv;
    }

    /// True if a key exists.
    bool has(const std::string& key) const
    {
        return _kv.find(key) != _kv.end();
    }

    /// Typed getter with default.
    template<class T>
    T get(const std::string& key, const T& def) const
    {
        auto it = _kv.find(key);
        if (it == _kv.end()) {
            return def;
        }
        return convert<T>(it->second);
    }

    /// Typed getter, throws if missing.
    template<class T>
    T require(const std::string& key) const
    {
        auto it = _kv.find(key);
        if (it == _kv.end()) {
            throw std::runtime_error("Missing required key: " + key);
        }
        return convert<T>(it->second);
    }

private:
    Scope        _scope;
    std::string  _name;
    Map          _kv;

    Keyword(Scope scope, std::string name, Map map)
        : _scope(std::move(scope)), _name(std::move(name)), _kv(std::move(map)) {}

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
