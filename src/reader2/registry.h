#pragma once
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include "command_spec.h"
#include "types.h"

namespace fem::reader2 {

    /**
     * \brief Global registry that stores command specifications by scope.
     */
    class Registry {
    public:
        using Configure = std::function<void(CommandSpec&)>;

        /// \brief Access the singleton registry instance.
        static Registry& instance();

        /// \brief Register a command handler with the registry.
        void register_command(Scope scope, std::string name, Configure cfg);
        /// \brief Locate a command specification for the given scope/name pair.
        const CommandSpec* find(const Scope& scope, const std::string& name) const;

        /// \brief Determine if a scope has registered children.
        bool scope_has_children(const Scope& scope) const;
        /// \brief Determine if a scope allows a particular child command.
        bool scope_allows(const Scope& scope, const std::string& child_name) const;
        /// \brief Collect all child command specifications for a scope.
        std::vector<const CommandSpec*> children_of(const Scope& scope_name) const;

        /// \brief Render documentation for a single command.
        std::string document(const Scope& scope, const std::string& name) const;
        /// \brief Render documentation for all commands.
        std::string document_all() const;
        /// \brief Print a tree representation of the registered scopes.
        void print_tree(std::ostream& os = std::cout) const;

    private:
        std::unordered_map<std::string, CommandSpec>            _map;
        std::unordered_map<Scope, std::unordered_set<std::string>> _children;
        static std::string key(const Scope& scope, const std::string& name);
    };

} // namespace fem::reader2
