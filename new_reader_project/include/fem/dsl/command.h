/**
 * @file command.h
 * @brief Declares `Command`, the top-level specification for a DSL keyword.
 *
 * A `Command` represents a registered keyword (e.g. `"ELASTIC"`, `"NUMEIGENVALUES"`).
 * Each command may:
 *
 *  - define an optional **admission condition** (`allow_if(...)`) that decides whether
 *    the command is allowed under the current parent scope (evaluated against
 *    `ParentInfo` and the command's own `Keys`);
 *  - provide one or more **variants** (`variant(...)`). The engine selects the first
 *    variant whose condition evaluates to true (or that has no condition) and then
 *    executes its segments.
 *
 * Typical usage:
 * @code
 * reg.command("ELASTIC", [](Command& c){
 *     c.allow_if( parent_is("MATERIAL") );
 *     c.variant( Variant::make()
 *         .segment( Segment::make()
 *             .range(LineRange{}.min(1).max(1))
 *             .pattern(Pattern::make().fixed<double,1>().fixed<double,1>())
 *             .bind([](double E, double nu){ ... })
 *         )
 *     );
 * });
 * @endcode
 *
 * Semantics:
 *  - If no admission condition is set (`allow_if` never called), `admit_` defaults to
 *    `Condition{}` which is the `Always` condition (i.e., admissible under any parent).
 *  - Variants are checked in the order they were added; the first matching one is used.
 *
 * @see variant.h
 * @see condition.h
 * @date 12.10.2025
 */

#pragma once
#include <string>
#include <utility>
#include <vector>
#include "condition.h"
#include "variant.h"

namespace fem {
namespace dsl {

/**
 * @class Command
 * @brief Top-level specification for a registered DSL keyword.
 */
struct Command {
    /** Command name as it appears in the input (normalized, e.g. `"ELASTIC"`). */
    std::string name_;

    /** Admission rule for this command; defaults to `Always`. */
    Condition admit_ = Condition{};

    /** Ordered list of variants; the engine picks the first matching one. */
    std::vector<Variant> variants_;

    /** Optional short description used by documentation printers. */
    std::string doc_;

    /**
     * @brief Constructs a command with the given name.
     */
    explicit Command(std::string n) : name_(std::move(n)) {}

    /**
     * @brief Sets the admission condition for this command.
     *
     * When set, the engine evaluates the condition against the current parent context
     * and the command's own keys. If the condition returns false, the engine climbs
     * the scope upwards until an ancestor admits the command or the scope is exhausted.
     *
     * @param c Condition tree to store by value (use helpers from `condition.h`).
     * @return Reference to `*this` for fluent chaining.
     */
    Command& allow_if(Condition c) {
        admit_ = std::move(c);
        return *this;
    }

    /**
     * @brief Appends a variant definition to this command.
     *
     * Variants are evaluated in insertion order; the first variant whose condition
     * holds (or has no condition) is selected and executed.
     *
     * @param v Variant to add (moved in).
     * @return Reference to `*this` for fluent chaining.
     */
    Command& variant(Variant v) {
        variants_.push_back(std::move(v));
        return *this;
    }

    /**
    * @brief Sets a human-readable description for this command (shown in help output).
    * @param d Short description text.
    * @return Reference to `*this` for fluent chaining.
    */
    Command& doc(std::string d) {
        doc_ = std::move(d);
        return *this;
    }
};

} // namespace dsl
} // namespace fem
