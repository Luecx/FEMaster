#pragma once
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
 * \brief Configuration container describing a registered command.
 */
struct CommandSpec {
    // Identity
    Scope       scope;
    std::string name;

    // Optional scope opening
    std::optional<Scope> open_scope;

    // Keyword validation and docs
    std::optional<KeyRules> key_rules;

    // Callbacks
    std::function<void(Context&, const Keyword&)> on_enter_cb;
    std::function<void(Context&, const Keyword&)> on_keyword_cb;
    std::function<void(Context&)>                 on_exit_cb;

    // Command-level docs
    std::string doc_title;
    std::string doc_summary;
    std::string doc_notes;

    /// \brief Document the command for generated reference material.
    CommandSpec& doc(std::string title, std::string summary, std::string notes = {})
    {
        doc_title   = std::move(title);
        doc_summary = std::move(summary);
        doc_notes   = std::move(notes);
        return *this;
    }

    // Fluent plans
    std::vector<ConditionalPlan> plans;

    // Convenience: add one plan from pieces
    /// \brief Add an execution plan for the command.
    CommandSpec& plan(Condition cnd,
                      std::initializer_list<Segment> segs,
                      std::string display_name = {})
    {
        ConditionalPlan p = ConditionalPlan::make()
                                .when(std::move(cnd))
                                .segments(segs)
                                .name(std::move(display_name));
        plans.push_back(std::move(p));
        return *this;
    }

    // Convenience setters
    /// \brief Declare that this command opens a child scope.
    CommandSpec& opens_scope(Scope s)
    {
        open_scope = std::move(s);
        return *this;
    }

    /// \brief Configure keyword validation rules.
    CommandSpec& keys(KeyRules k)
    {
        key_rules = std::move(k);
        return *this;
    }

    /// \brief Register an on-enter callback.
    CommandSpec& on_enter(std::function<void(Context&, const Keyword&)> f)
    {
        on_enter_cb = std::move(f);
        return *this;
    }

    /// \brief Register an on-keyword callback.
    CommandSpec& on_keyword(std::function<void(Context&, const Keyword&)> f)
    {
        on_keyword_cb = std::move(f);
        return *this;
    }

    /// \brief Register an on-exit callback.
    CommandSpec& on_exit(std::function<void(Context&)> f)
    {
        on_exit_cb = std::move(f);
        return *this;
    }
};

} // namespace fem::reader2
