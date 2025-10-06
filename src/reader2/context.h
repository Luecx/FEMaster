#pragma once
#include <vector>
#include <stdexcept>
#include <numeric>
#include "types.h"

namespace fem::reader2 {

    /**
     * @brief Tracks the current parser scope hierarchy (e.g. ROOT → MATERIAL → ELASTIC).
     */
    struct Context {
        /// Construct with a configurable root scope (defaults to "ROOT").
        explicit Context(Scope root_name = "ROOT")
            : scope_stack{std::move(root_name)} {}

        /// Current active scope (top of stack).
        const Scope& current_scope() const { return scope_stack.back(); }

        /// Root scope name (bottom of stack).
        const Scope& root() const { return scope_stack.front(); }

        /// True if currently at the root scope.
        bool at_root() const { return scope_stack.size() == 1; }

        /// Push a new scope.
        void push_scope(Scope s) { scope_stack.push_back(std::move(s)); }

        /// Pop the current scope, but never below the root.
        void pop_scope() {
            if (at_root())
                throw std::runtime_error("Context::pop_scope(): cannot pop the root scope");
            scope_stack.pop_back();
        }

        /// Returns the full scope path, like "ROOT/MATERIAL/ELASTIC".
        std::string path(char sep = '/') const {
            std::string out;
            for (size_t i = 0; i < scope_stack.size(); ++i) {
                if (i) out.push_back(sep);
                out += scope_stack[i];
            }
            return out;
        }

        std::vector<Scope> scope_stack;
    };

} // namespace fem::reader2
