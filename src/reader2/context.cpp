/**
* @file context.cpp
 * @brief Implements the Context struct for tracking parser scope hierarchy.
 *
 * @date 06.10.2025
 * @version 1.0
 */

#include "context.h"
#include <stdexcept>

namespace fem::reader2 {

    Context::Context(Scope root_name)
        : _scopeStack{std::move(root_name)} {}

    const Scope& Context::current_scope() const {
        return _scopeStack.back();
    }

    const Scope& Context::root() const {
        return _scopeStack.front();
    }

    bool Context::at_root() const {
        return _scopeStack.size() == 1;
    }

    size_t Context::depth() const {
        return _scopeStack.size();
    }

    void Context::push_scope(Scope s) {
        _scopeStack.push_back(std::move(s));
    }

    void Context::pop_scope() {
        if (at_root()) {
            throw std::runtime_error("Context::pop_scope(): cannot pop the root scope");
        }
        _scopeStack.pop_back();
    }

    std::string Context::path(char sep) const {
        std::string out;
        out.reserve(32); // small optimization; grows as needed
        for (size_t i = 0; i < _scopeStack.size(); ++i) {
            if (i) out.push_back(sep);
            out += _scopeStack[i];
        }
        return out;
    }

} // namespace fem::reader2
