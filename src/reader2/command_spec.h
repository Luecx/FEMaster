// parser/command_spec.h
#pragma once
/**
 * @file command_spec.h
 * @brief Declarative command specification: schemas, callbacks, and scope behavior.
 */

#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "context.h"
#include "schema.h"
#include "key_rules.h"
#include "types.h"

namespace fem::reader2 {

// Forward declaration to avoid heavy includes in headers using CommandSpec.
class Keyword;

/**
 * @brief Specification for a single command within a parent scope.
 *
 * A CommandSpec describes:
 *  - Which scope it belongs to (parent scope)
 *  - Its name
 *  - Whether it opens a child scope
 *  - Optional key-validation rules
 *  - Callbacks on enter/keyword/exit
 *  - Data-line handling mode: simple lines, conditional lines, or token collection
 */
struct CommandSpec {
    // Identity in registry: (scope, name)
    Scope       scope;
    std::string name;

    // If the command has children, a scope opens automatically (default: use `name`)
    std::optional<Scope> open_scope;

    // Keyword key validation and defaults/enums
    std::optional<KeyRules> key_rules;

    // Optional callbacks
    std::function<void(Context&, const Keyword&)> on_enter_cb;   // before data processing
    std::function<void(Context&, const Keyword&)> on_keyword_cb; // immediate
    std::function<void(Context&)>                 on_exit_cb;    // when this scope auto-closes

    // --- Documentation fields (command-level) ---
    std::string doc_title;   ///< e.g. "Support"
    std::string doc_summary; ///< short paragraph describing command intent
    std::string doc_notes;   ///< extra free-form notes (constraints, tips, pitfalls)

    CommandSpec& doc(std::string title, std::string summary, std::string notes = {}){
        doc_title   = std::move(title);
        doc_summary = std::move(summary);
        doc_notes   = std::move(notes);
        return *this;
    }

    // --- Variant A: simple fixed- or variable-line schemas (what you already have) ---
    bool has_lines = false;
    LineRange range{0,0};
    Schema schema;

    // --- Variant B: conditional schemas selected from the keyword args ---
    struct ConditionalLines {
        std::function<bool(const Keyword&)> when;  // guard (now sees Keyword)
        LineRange range;
        Schema schema;

        // doc block for this alternative
        std::string title;       // label for the branch, e.g. "TYPE=NODE"
        std::string description; // human explanation
    };
    std::vector<ConditionalLines> conditional_lines; // first matching guard wins

    // doc helper for a conditional branch
    CommandSpec& lines_if(std::function<bool(const Keyword&)> when, LineRange r, const Schema& s,
                          std::string title = {}, std::string desc = {}){
        conditional_lines.push_back({std::move(when), r, s, std::move(title), std::move(desc)});
        return *this;
    }

    // --- Variant C: token accumulation across multiple data lines ---
    struct TokenCollect {
        std::function<size_t(const Keyword&)> required;
        std::function<void(Context&, const std::vector<std::string>&, const LineMeta&)> handle;

        // docs
        std::string title;
        std::string detail;  // include semantics, order, examples
        std::vector<std::string> token_labels;  // optional labels (e.g. ["EID","N1..N20"])
    };
    std::optional<TokenCollect> token_collect;

    CommandSpec& collect_tokens(std::function<size_t(const Keyword&)> required,
                                std::function<void(Context&, const std::vector<std::string>&, const LineMeta&)> handler,
                                std::string title = {}, std::string detail = {},
                                std::initializer_list<std::string> labels = {}){
        token_collect = TokenCollect{ std::move(required), std::move(handler),
                                      std::move(title), std::move(detail), {labels.begin(), labels.end()} };
        return *this;
    }

    // ---- Fluent API (rest) ----
    CommandSpec& opens_scope(Scope s){ open_scope = std::move(s); return *this; }
    CommandSpec& keys(KeyRules k){ key_rules = std::move(k); return *this; }
    CommandSpec& on_enter(std::function<void(Context&, const Keyword&)> f){ on_enter_cb = std::move(f); return *this; }
    CommandSpec& on_keyword(std::function<void(Context&, const Keyword&)> f){ on_keyword_cb = std::move(f); return *this; }
    CommandSpec& on_exit(std::function<void(Context&)> f){ on_exit_cb = std::move(f); return *this; }
    CommandSpec& lines(LineRange r, const Schema& s){ has_lines=true; range=r; schema=s; return *this; }
};

} // namespace fem::reader2
