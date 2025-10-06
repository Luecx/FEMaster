// parser/reader.cpp
/**
 * @file reader.cpp
 * @brief Streaming reader that walks a File, resolves scopes, and dispatches to CommandSpec.
 */
#include "reader.h"
#include "registry.h"
#include "types.h"
#include "keyword.h"

#include <stdexcept>
#include <sstream>
#include <iostream>

namespace fem::reader2 {

//------------------------------------------------------------------------------
// Helper: next_effective / peek_effective
//------------------------------------------------------------------------------
Line& Reader::next_effective(File& f) {
    if (look_) {
        scratch_ = *look_;
        look_.reset();
        return scratch_;
    }
    return f.next_effective();
}

Line& Reader::peek_effective(File& f) {
    if (!look_) look_ = f.next_effective();
    return *look_;
}

//------------------------------------------------------------------------------
// ensure_scope_accepts_or_bubble()
// - Checks whether current scope (or its ancestors) accept a command.
// - If not accepted, automatically closes scopes upward, calling on_exit.
// - Returns index in scope_stack of accepting scope.
//------------------------------------------------------------------------------
size_t Reader::ensure_scope_accepts_or_bubble(const std::string& name) {
    while (true) {
        const Scope& scope = cx_.current_scope();

        // 1) Does current scope handle this command?
        if (Registry::instance().find(scope, name)) {
            return cx_.scope_stack.size() - 1;
        }

        // 2) If we’re at ROOT and still no match → error
        if (scope == "ROOT") {
            throw std::runtime_error(
                "Unexpected command '" + name +
                "' at ROOT (no matching handler registered).");
        }

        // 3) Otherwise: close this scope properly and bubble up
        if (!blocks_.empty() && blocks_.back().scope_name == scope) {
            auto blk = blocks_.back();
            blocks_.pop_back();
            if (blk.on_exit_cb)
                blk.on_exit_cb(cx_);
        }

        cx_.pop_scope();  // move one level up and retry
    }
}

//------------------------------------------------------------------------------
// run()
//------------------------------------------------------------------------------
void Reader::run(const std::string& filepath) {
    File file(filepath);

    while (true) {
        Line& L = peek_effective(file);
        if (L.type() == LineType::EOF_MARK) break;

        if (L.type() != LineType::KEYWORD) {
            look_.reset();
            continue;
        }

        std::string name = L.command();

        // Automatically close any scopes that do not accept this command.
        ensure_scope_accepts_or_bubble(name);

        const Scope& scope = cx_.current_scope();
        const CommandSpec* spec = Registry::instance().find(scope, name);
        if (!spec) {
            look_.reset();
            continue;
        }

        // Build Keyword (parsing + KeyRules application centralized)
        Keyword kw = Keyword::from_line(L, scope, spec->key_rules);

        // Consume this keyword line
        scratch_ = *look_;
        look_.reset();

        // Callbacks for entering
        if (spec->on_enter_cb)
            spec->on_enter_cb(cx_, kw);
        if (spec->on_keyword_cb)
            spec->on_keyword_cb(cx_, kw);

        //===========================================================
        // MODE C: token accumulation (e.g. collect N numbers)
        //===========================================================
        if (spec->token_collect) {
            const size_t need = spec->token_collect->required(kw);
            std::vector<std::string> tokens;
            tokens.reserve(need);

            LineMeta last_meta{file.path(), file.line_number()};

            while (tokens.size() < need) {
                Line& P = peek_effective(file);
                if (P.type() == LineType::EOF_MARK) break;
                if (P.type() != LineType::DATA) break;

                last_meta = {file.path(), file.line_number()};
                const auto& vals = P.values();
                tokens.insert(tokens.end(), vals.begin(), vals.end());

                scratch_ = *look_;
                look_.reset();
            }

            if (tokens.size() < need) {
                throw std::runtime_error("Insufficient values for '" + name +
                    "': required " + std::to_string(need) +
                    ", got " + std::to_string(tokens.size()));
            }
            if (tokens.size() > need) {
                throw std::runtime_error("Too many values for '" + name +
                    "': required " + std::to_string(need) +
                    ", got " + std::to_string(tokens.size()));
            }

            spec->token_collect->handle(cx_, tokens, last_meta);
        }

        //===========================================================
        // MODE B: conditional schemas
        //===========================================================
        else if (!spec->conditional_lines.empty()) {
            const CommandSpec::ConditionalLines* chosen = nullptr;
            for (auto& v : spec->conditional_lines) {
                if (v.when(kw)) {
                    chosen = &v;
                    break;
                }
            }
            if (!chosen) {
                throw std::runtime_error("No conditional schema matched for '" + name + "'");
            }

            size_t consumed = 0;
            while (consumed < chosen->range.max_lines) {
                Line& P = peek_effective(file);
                if (P.type() == LineType::EOF_MARK) break;
                if (P.type() == LineType::KEYWORD) break;
                if (P.type() != LineType::DATA) {
                    look_.reset();
                    break;
                }

                bool matched = false;
                for (auto& alt : chosen->schema.alternatives) {
                    LineMeta meta{file.path(), file.line_number()};
                    if (alt(cx_, P.values(), meta)) {
                        matched = true;
                        break;
                    }
                }
                if (!matched) {
                    throw std::runtime_error("No schema matched at " + file.path() +
                        ":" + std::to_string(file.line_number()));
                }

                scratch_ = *look_;
                look_.reset();
                ++consumed;
            }
            if (consumed < chosen->range.min_lines) {
                throw std::runtime_error("Too few data lines for '" + name +
                    "' (min=" + std::to_string(chosen->range.min_lines) + ")");
            }
        }

        //===========================================================
        // MODE A: simple schema
        //===========================================================
        else if (spec->has_lines) {
            size_t consumed = 0;
            while (consumed < spec->range.max_lines) {
                Line& P = peek_effective(file);
                if (P.type() == LineType::EOF_MARK) break;
                if (P.type() == LineType::KEYWORD) break;
                if (P.type() != LineType::DATA) {
                    look_.reset();
                    break;
                }

                bool matched = false;
                for (auto& alt : spec->schema.alternatives) {
                    LineMeta meta{file.path(), file.line_number()};
                    if (alt(cx_, P.values(), meta)) {
                        matched = true;
                        break;
                    }
                }
                if (!matched) {
                    throw std::runtime_error("No schema matched at " + file.path() +
                        ":" + std::to_string(file.line_number()));
                }

                scratch_ = *look_;
                look_.reset();
                ++consumed;
            }
            if (consumed < spec->range.min_lines) {
                throw std::runtime_error("Too few data lines for '" + name +
                    "' (min=" + std::to_string(spec->range.min_lines) + ")");
            }
        }

        //===========================================================
        // Scope management
        //===========================================================
        const Scope open_as = spec->open_scope ? *spec->open_scope : spec->name;

        bool explicit_scope = static_cast<bool>(spec->open_scope);
        bool has_children   = Registry::instance().scope_has_children(open_as);

        // Always open explicitly declared scopes (like MATERIAL)
        if (explicit_scope) {
            cx_.push_scope(open_as);
            blocks_.push_back(ActiveBlock{open_as, spec->on_exit_cb});
        }
        // Automatically open scopes that have known children
        else if (has_children) {
            cx_.push_scope(open_as);
            blocks_.push_back(ActiveBlock{open_as, spec->on_exit_cb});
        }
    }

    // --- Ensure all blocks are closed at EOF ---
    while (!blocks_.empty()) {
        auto blk = blocks_.back();
        blocks_.pop_back();
        if (blk.on_exit_cb)
            blk.on_exit_cb(cx_);
        cx_.scope_stack.pop_back();
    }
}

} // namespace fem::reader2
