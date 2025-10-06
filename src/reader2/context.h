#pragma once
/**
 * @file context.h
 * @brief Declares the Context struct that tracks the current parser scope hierarchy
 *        (e.g., ROOT → MATERIAL → ELASTIC).
 *
 * @date 06.10.2025
 * @version 1.0
 */

#include <string>
#include <vector>
#include "types.h"

namespace fem::reader2 {

/**
 * @brief Tracks the current parser scope hierarchy.
 */
struct Context {
    /** @brief Construct with a configurable root scope (defaults to "ROOT"). */
    explicit Context(Scope root_name = "ROOT");

    /** @brief Current active scope (top of stack). */
    const Scope& current_scope() const;

    /** @brief Root scope name (bottom of stack). */
    const Scope& root() const;

    /** @brief True if currently at the root scope. */
    bool at_root() const;

    /** @brief Current depth of the scope stack. */
    size_t depth() const;

    /** @brief Push a new scope. */
    void push_scope(Scope s);

    /**
     * @brief Pop the current scope, but never below the root.
     * @throws std::runtime_error if called at root.
     */
    void pop_scope();

    /**
     * @brief Returns the full scope path, like "ROOT/MATERIAL/ELASTIC".
     * @param sep Path separator character (default '/').
     */
    std::string path(char sep = '/') const;

    std::vector<Scope> _scopeStack; ///< Internal stack, front = root, back = current
};

} // namespace fem::reader2
