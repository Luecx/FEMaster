// parser/registry.h
#pragma once
/**
 * @file registry.h
 * @brief Global command registry and documentation helpers.
 */

#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>

#include "command_spec.h"
#include "types.h"

namespace fem::reader2 {

    class Registry {
    public:
        using Configure = std::function<void(CommandSpec&)>;

        static Registry& instance();

        void register_command(Scope scope, std::string name, Configure cfg);
        const CommandSpec* find(const Scope& scope, const std::string& name) const;

        bool scope_has_children(const Scope& scope) const;
        bool scope_allows(const Scope& scope, const std::string& child_name) const;
        std::vector<const CommandSpec*> children_of(const Scope& scope_name) const;

        std::string document(const Scope& scope, const std::string& name) const;
        std::string document_all() const;

        /// Print all commands as a scope tree to std::ostream (default = std::cout)
        void print_tree(std::ostream& os = std::cout) const;

    private:
        std::unordered_map<std::string, CommandSpec> map_;
        std::unordered_map<Scope, std::unordered_set<std::string>> children_;

        static std::string key(const Scope& scope, const std::string& name);
    };

    struct CommandRegistrar {
        CommandRegistrar(const char* scope, const char* name, Registry::Configure cfg);
    };

} // namespace fem::reader2
