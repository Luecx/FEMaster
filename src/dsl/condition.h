#pragma once
/**
 * @file condition.h
 * @brief Declarative condition tree to admit commands and select variants in the DSL.
 *
 * This header is self-contained: it defines the Condition type, its factory helpers
 * as static methods, the evaluation logic, and (optional) free-function forwarders
 * for backwards compatibility.
 *
 * ## Overview
 * A Condition is a small AST of predicates evaluated against:
 *   - the current parent-scope context (`ParentInfo`), and
 *   - the command's own keyword-line keys (`Keys`).
 *
 * Typical usage with static factories on Condition:
 * @code
 * using fem::dsl::Condition;
 *
 * c.allow_if( Condition::parent_is("MATERIAL") );
 *
 * c.allow_if( Condition::all_of({
 *   Condition::parent_is("LOADCASE"),
 *   Condition::parent_key_equals("TYPE", {"LINEARBUCKLING","EIGENFREQ"})
 * }) );
 *
 * v.when( Condition::key_equals("MODE", {"A","B"}) );
 * @endcode
 *
 * If you have existing code using free helpers like `parent_is(...)` or `all_of(...)`,
 * this header also provides inline forwarders to the new `Condition::...` APIs.
 *
 * ## Kinds
 *  - Structural:   All / Any / Not / Always
 *  - Parent-based: ParentCommandIs / ParentKeyEquals / ParentHasKey
 *  - Self-based:   KeyPresent / KeyEquals
 *
 * ## Evaluation
 * `cond.eval(parent, self_keys)` returns true/false by recursively evaluating the tree.
 *
 * @date 12.10.2025
 */

#include <initializer_list>
#include <string>
#include <vector>
#include <cstddef>
#include "keys.h"

namespace fem { namespace dsl {

/**
 * @struct ParentInfo
 * @brief Lightweight view of the current parent in the scope stack.
 *
 * - `command`: name of the parent command (normalized, e.g. "MATERIAL")
 * - `keys`   : its keyword-line keys (as parsed by `Keys::from_keyword_line`)
 */
struct ParentInfo {
    std::string command;
    Keys        keys;
};

/**
 * @class Condition
 * @brief Tree-structured predicate evaluated against `(parent, self_keys)`.
 *
 * Build conditions via static factory methods:
 *  - Structure:
 *      - `Condition::always_true()`
 *      - `Condition::all_of({...})`
 *      - `Condition::any_of({...})`
 *      - `Condition::negate(child)`
 *  - Parent context:
 *      - `Condition::parent_is("MATERIAL")`
 *      - `Condition::parent_key_equals("TYPE", {"EIGENFREQ","LINEARBUCKLING"})`
 *      - `Condition::parent_has_key("NAME")`
 *  - Self keys (keys on the *current* keyword line):
 *      - `Condition::key_present("MODE")`
 *      - `Condition::key_equals("MODE", {"A","B"})`
 *
 * Example:
 * @code
 * Condition c = Condition::all_of({
 *   Condition::parent_is("LOADCASE"),
 *   Condition::parent_key_equals("TYPE", {"EIGENFREQ","LINEARBUCKLING"}),
 *   Condition::negate( Condition::key_equals("DEBUG","ON") )
 * });
 * @endcode
 */
struct Condition {
    /// Node kinds for the condition tree.
    enum Kind {
        Always,            ///< Always true
        All,               ///< Conjunction over `children`
        Any,               ///< Disjunction over `children`
        Not,               ///< Logical negation of first/only child
        ParentCommandIs,   ///< `parent.command` is in `values`
        ParentKeyEquals,   ///< `parent.keys.raw(key)` equals any of `values`
        ParentHasKey,      ///< `parent.keys.has(key)`
        KeyPresent,        ///< `self_keys.has(key)`
        KeyEquals          ///< `self_keys.raw(key)` equals any of `values`
    };

    // --- Data ---
    Kind kind = Always;                      ///< Node kind (defaults to `Always`).
    std::string key;                         ///< Optional key (for key-based kinds).
    std::vector<std::string> values;         ///< Optional value set.
    std::vector<Condition> children;         ///< Child nodes (for All/Any/Not).

    /// Default constructor → `Always`.
    Condition() = default;

    // ---------------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------------

    /**
     * @brief Evaluates this condition against the given context.
     * @param parent    Context of the current parent in the scope stack.
     * @param self_keys Keys from the current command's keyword line.
     * @return `true` if the condition holds, otherwise `false`.
     */
    bool eval(const ParentInfo& parent, const Keys& self_keys) const {
        switch (kind) {
            case Always:
                return true;

            case All:
                for (const auto& ch : children)
                    if (!ch.eval(parent, self_keys)) return false;
                return true;

            case Any:
                for (const auto& ch : children)
                    if (ch.eval(parent, self_keys)) return true;
                return false;

            case Not:
                if (children.empty()) return true;
                return !children.front().eval(parent, self_keys);

            case ParentCommandIs:
                for (const auto& v : values)
                    if (v == parent.command) return true;
                return false;

            case ParentKeyEquals: {
                const std::string& raw = parent.keys.raw(key);
                if (raw.empty()) return false;
                for (const auto& v : values)
                    if (v == raw) return true;
                return false;
            }

            case ParentHasKey:
                return parent.keys.has(key);

            case KeyPresent:
                return self_keys.has(key);

            case KeyEquals: {
                const std::string& raw = self_keys.raw(key);
                if (raw.empty()) return false;
                for (const auto& v : values)
                    if (v == raw) return true;
                return false;
            }
        }
        return false; // unreachable, but keeps some compilers happy
    }

    // ---------------------------------------------------------------------
    // Static factory helpers — STRUCTURE
    // ---------------------------------------------------------------------

    /**
     * @brief Always-true condition.
     */
    static Condition always_true() {
        return Condition{};
    }

    /**
     * @brief Conjunction (logical AND) over a list of child conditions.
     */
    static Condition all_of(std::initializer_list<Condition> conds) {
        Condition c;
        c.kind = All;
        c.children.assign(conds.begin(), conds.end());
        return c;
    }

    /**
     * @brief Disjunction (logical OR) over a list of child conditions.
     */
    static Condition any_of(std::initializer_list<Condition> conds) {
        Condition c;
        c.kind = Any;
        c.children.assign(conds.begin(), conds.end());
        return c;
    }

    /**
     * @brief Logical negation of a child condition.
     */
    static Condition negate(const Condition& cond) {
        Condition c;
        c.kind = Not;
        c.children.clear();
        c.children.push_back(cond);
        return c;
    }

    // ---------------------------------------------------------------------
    // Static factory helpers — PARENT CONTEXT
    // ---------------------------------------------------------------------

    /**
     * @brief Matches the parent command name against one or more expected names.
     * @code
     * Condition::parent_is("MATERIAL")
     * Condition::parent_is({"MATERIAL","SECTION"})
     * @endcode
     */
    static Condition parent_is(std::initializer_list<std::string> names) {
        Condition c;
        c.kind = ParentCommandIs;
        c.values.assign(names.begin(), names.end());
        return c;
    }

    /// @overload Convenience overload for single name.
    static Condition parent_is(const std::string& name) {
        return parent_is({name});
    }

    /**
     * @brief Checks that the parent's key `k` equals any of `vals`.
     * The comparison uses the raw (already-normalized) string value.
     */
    static Condition parent_key_equals(const std::string& k,
                                       std::initializer_list<std::string> vals) {
        Condition c;
        c.kind = ParentKeyEquals;
        c.key = k;
        c.values.assign(vals.begin(), vals.end());
        return c;
    }

    /// @overload Single-value convenience.
    static Condition parent_key_equals(const std::string& k, const std::string& value) {
        return parent_key_equals(k, {value});
    }

    /**
     * @brief Requires that the parent has a given key, regardless of its value.
     */
    static Condition parent_has_key(const std::string& k) {
        Condition c;
        c.kind = ParentHasKey;
        c.key = k;
        return c;
    }

    // ---------------------------------------------------------------------
    // Static factory helpers — SELF KEYS (current keyword line)
    // ---------------------------------------------------------------------

    /**
     * @brief Requires that the current command line contains the key `k`.
     */
    static Condition key_present(const std::string& k) {
        Condition c;
        c.kind = KeyPresent;
        c.key = k;
        return c;
    }

    /**
     * @brief Requires that the current command line's key `k` equals any of `vals`.
     */
    static Condition key_equals(const std::string& k,
                                std::initializer_list<std::string> vals) {
        Condition c;
        c.kind = KeyEquals;
        c.key = k;
        c.values.assign(vals.begin(), vals.end());
        return c;
    }

    /// @overload Single-value convenience.
    static Condition key_equals(const std::string& k, const std::string& value) {
        return key_equals(k, {value});
    }
};

}} // namespace fem::dsl