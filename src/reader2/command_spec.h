#pragma once
/**
 * @file command_spec.h
 * @brief Declares the CommandSpec struct, describing a registered command
 *        within the FEMaster parser framework.
 *
 * @details
 * A CommandSpec represents a single command keyword and defines its behavior,
 * documentation, callbacks, and conditional execution plans. It can also
 * open new parsing scopes and define key validation rules.
 *
 * @date 06.10.2025
 * @version 1.0
 * @see conditional_plan.h
 */

#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "context.h"
#include "pattern.h"
#include "segment.h"
#include "conditional_plan.h"
#include "line_range.h"
#include "key_rules.h"
#include "types.h"
#include "condition.h"

namespace fem::reader2 {

class Keyword;

/**
 * @brief Configuration container describing a registered parser command.
 */
struct CommandSpec {
    // Identity
    Scope       scope;              ///< Scope in which the command is valid
    std::string name;               ///< Command keyword name
    std::optional<Scope> open_scope; ///< Optional child scope

    // Keyword validation
    std::optional<KeyRules> key_rules; ///< Validation rules for keyword arguments

    // Callbacks
    std::function<void(Context&, const Keyword&)> on_enter_cb;
    std::function<void(Context&, const Keyword&)> on_keyword_cb;
    std::function<void(Context&)>                 on_exit_cb;

    // Documentation
    std::string doc_summary;
    std::string doc_notes;

    // Execution plans
    std::vector<ConditionalPlan> plans;

    /** @brief Document the command for generated reference material. */
    CommandSpec& doc(std::string summary, std::string notes = {});

    /** @brief Add an execution plan for the command. */
    CommandSpec& plan(Condition cnd,
                      std::initializer_list<Segment> segs,
                      std::string display_name = {});

    /** @brief Add an execution plan consisting of a single segment. */
    CommandSpec& plan(Condition cnd,
                      Segment segment,
                      std::string display_name = {});

    /** @brief Declare that this command opens a child scope. */
    CommandSpec& opens_scope(Scope s);

    /** @brief Configure keyword validation rules. */
    CommandSpec& keys(KeyRules k);

    /** @brief Register an on-enter callback. */
    CommandSpec& on_enter(std::function<void(Context&, const Keyword&)> f);

    /** @brief Register an on-keyword callback. */
    CommandSpec& on_keyword(std::function<void(Context&, const Keyword&)> f);

    /** @brief Register an on-exit callback. */
    CommandSpec& on_exit(std::function<void(Context&)> f);
};

} // namespace fem::reader2
