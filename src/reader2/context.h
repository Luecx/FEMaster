#pragma once
#include <stdexcept>
#include <vector>

#include "types.h"

namespace fem::reader2 {

    /**
     * @brief Tracks the current parser scope hierarchy (e.g. ROOT → MATERIAL → ELASTIC).
     */
    struct Context {
        /// Construct with a configurable root scope (defaults to "ROOT").
        explicit Context(Scope root_name = "ROOT")
            : _scopeStack{std::move(root_name)} {}

        /// Current active scope (top of stack).
        const Scope& current_scope() const
        {
            return _scopeStack.back();
        }

        /// Root scope name (bottom of stack).
        const Scope& root() const
        {
            return _scopeStack.front();
        }

        /// True if currently at the root scope.
        bool at_root() const
        {
            return _scopeStack.size() == 1;
        }

        /// Current depth of the scope stack.
        size_t depth() const
        {
            return _scopeStack.size();
        }

        /// Push a new scope.
        void push_scope(Scope s)
        {
            _scopeStack.push_back(std::move(s));
        }

        /// Pop the current scope, but never below the root.
        void pop_scope()
        {
            if (at_root())
                throw std::runtime_error("Context::pop_scope(): cannot pop the root scope");
            _scopeStack.pop_back();
        }

        /// Returns the full scope path, like "ROOT/MATERIAL/ELASTIC".
        std::string path(char sep = '/') const
        {
            std::string out;
            for (size_t i = 0; i < _scopeStack.size(); ++i) {
                if (i) out.push_back(sep);
                out += _scopeStack[i];
            }
            return out;
        }

        std::vector<Scope> _scopeStack;
    };

} // namespace fem::reader2
