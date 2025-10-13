#pragma once
/**
 * @file scope_condition.h
 * @brief Boolescher Baum für Scoping-Prüfungen: Orts-Atome (At/Within) + optionale KeyCondition,
 *        kombinierbar über And/Or/Not.
 */

#include "context.h"
#include "key_condition.h"
#include "keyword.h"
#include "types.h"

#include <initializer_list>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

namespace fem::reader2::condition {

/**
 * @brief Scoping-Bedingung (At/Within + KeyCondition) mit Bool-Kombinatoren.
 */
class ScopeCondition {
    public:
    enum class Op { Atom, And, Or, Not };
    enum class Location { At, Within };

    // Atom-Fabriken
    static ScopeCondition At(std::initializer_list<fem::reader2::Scope> scopes);
    static ScopeCondition Within(std::initializer_list<fem::reader2::Scope> scopes);

    // Optional: KeyCondition an Atom hängen
    ScopeCondition& Where(KeyCondition expr);

    // Bool-Kombinatoren
    static ScopeCondition And(std::initializer_list<ScopeCondition> xs);
    static ScopeCondition Or (std::initializer_list<ScopeCondition> xs);
    static ScopeCondition Not(ScopeCondition x);

    bool satisfied_by(const fem::reader2::Context& ctx, bool nearest_first = true) const;

    private:
    struct Atom {
        Location                             loc = Location::At;
        std::unordered_set<fem::reader2::Scope> scopes;
        std::optional<KeyCondition>          expr;
    };

    Op                      _op { Op::Atom };
    Atom                    _atom {};
    std::vector<ScopeCondition> _children;

    bool eval_atom(const fem::reader2::Context& ctx, bool nearest_first) const;
};

} // namespace fem::reader2::condition
