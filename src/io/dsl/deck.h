/**
 * @file deck.h
 * @brief Defines the parsed input-deck representation used by FEMaster readers.
 *
 * `Deck` stores the syntactic result of one complete DSL parse. Each keyword
 * occurrence becomes one `ParsedCommand` with normalized keyword arguments,
 * its selected parent scope and the already validated segment invocations that
 * belong to its data block. Data lines are therefore not represented as tree
 * nodes; they remain grouped below their owning command exactly as required by
 * the selected DSL variant.
 *
 * The hierarchy has one real `ROOT` command at node zero. Every parsed command,
 * including top-level commands, therefore has a valid parent node and all tree
 * access uses the same `children()` and `execute_children()` operations.
 * Parsing and semantic execution remain intentionally separated: `DeckParser`
 * determines scopes, variants and record layouts once, while the higher-level
 * FEMaster reader chooses the dependency order in which stored commands execute.
 *
 * The stored command and segment pointers refer to the `Registry` used to parse
 * the deck. That registry must therefore outlive the `Deck` and all semantic
 * processing of its commands.
 *
 * @see DeckParser
 * @see Command
 * @see Segment
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "../../core/logging.h"
#include "command.h"
#include "line.h"

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::io::dsl {

/**
 * @brief One validated invocation of a DSL segment callback.
 *
 * The parser stores the normalized token sequence after pattern completion has
 * succeeded. Conversion to the bound C++ argument types remains owned by the
 * existing `Segment::_invoke` callback and therefore happens only when the
 * semantic reader explicitly executes this invocation.
 */
struct ParsedInvocation {
    const Segment*            segment = nullptr;
    std::vector<std::string> tokens;
    SourceLocation            location;
};

/**
 * @brief One keyword occurrence in a parsed input deck.
 *
 * Commands form a compact tree through integer parent/child indices. A node owns
 * only command-level state and its logical data-record invocations; individual
 * physical data lines never become tree nodes. Node zero is the synthetic ROOT
 * command and is its own parent, so every stored parent id is always valid.
 * Closing keywords such as `ENDPART` and `ENDSTEP` remain ordinary children of
 * the scope they terminate so their callbacks stay available during processing.
 */
struct ParsedCommand {
    using NodeId = std::size_t;

    const Command* command = nullptr;
    Keys           keys;

    NodeId parent = 0;

    std::vector<NodeId>           children;
    std::vector<ParsedInvocation> invocations;

    SourceLocation location;
};

/**
 * @brief Immutable parsed deck with explicit semantic execution operations.
 *
 * `Deck` owns one real ROOT node followed by the ordered command occurrences
 * produced by `DeckParser`. Queries return stable integer node ids, while
 * `enter()`, `leave()` and `execute()` invoke callbacks already registered on
 * the corresponding DSL command. Semantic readers can therefore state their
 * dependency order directly without maintaining a parallel root collection or
 * sentinel parent id.
 */
class Deck {
public:
    using NodeId = ParsedCommand::NodeId;

    Deck() {
        ParsedCommand root_node;
        root_node.command = &root_command();
        root_node.parent  = root();
        nodes_.push_back(std::move(root_node));
    }

    static constexpr NodeId root() {
        return 0;
    }

    // Access parsed command occurrences and filtered direct-child ranges
    const ParsedCommand& node(NodeId id) const {
        return nodes_.at(id);
    }

    std::vector<NodeId> children(NodeId parent, const std::string& command = {}) const {
        const auto& source = nodes_.at(parent).children;

        std::vector<NodeId> result;
        for (const NodeId id : source) {
            if (command.empty() || nodes_[id].command->name_ == command) {
                result.push_back(id);
            }
        }
        return result;
    }

    // Execute one stored command without implicitly traversing its children
    void enter(NodeId id) const {
        const ParsedCommand& current = nodes_.at(id);
        logging::info("Processing command: *", current.command->name_);

        if (current.command->on_enter_) {
            try {
                current.command->on_enter_(parent_info(current), current.keys);
            } catch (const std::exception& e) {
                throw_at(current.location, e.what());
            }
        }

        for (const ParsedInvocation& invocation : current.invocations) {
            if (!invocation.segment || !invocation.segment->_invoke) continue;
            try {
                invocation.segment->_invoke(invocation.tokens);
            } catch (const std::exception& e) {
                throw_at(invocation.location, e.what());
            }
        }
    }

    void leave(NodeId id) const {
        const ParsedCommand& current = nodes_.at(id);
        if (!current.command->on_exit_) return;

        try {
            current.command->on_exit_(current.keys);
        } catch (const std::exception& e) {
            throw_at(current.location, e.what());
        }
    }

    void execute(NodeId id) const {
        enter(id);
        leave(id);
    }

    // Execute selected direct children in the requested semantic order
    void execute_children(NodeId parent, const std::string& command) const {
        for (const NodeId id : children(parent, command)) execute(id);
    }

    void execute_children(NodeId parent, std::initializer_list<const char*> commands) const {
        for (const char* command : commands) execute_children(parent, command);
    }

    void execute_children(NodeId parent) const {
        for (const NodeId id : children(parent)) execute(id);
    }

    // Execute matching child scopes while keeping their callbacks active around
    // either the selected child commands or, for an empty list, all direct children.
    void execute_children(NodeId parent,
                          const std::string& command,
                          std::initializer_list<const char*> child_commands) const {
        for (const NodeId id : children(parent, command)) {
            enter(id);
            if (child_commands.size() == 0) {
                execute_children(id);
            } else {
                execute_children(id, child_commands);
            }
            leave(id);
        }
    }

private:
    std::vector<ParsedCommand> nodes_;

    friend class DeckParser;

    static const Command& root_command() {
        static const Command command{"ROOT"};
        return command;
    }

    ParentInfo parent_info(const ParsedCommand& current) const {
        const ParsedCommand& parent = nodes_.at(current.parent);
        return ParentInfo{parent.command->name_, parent.keys};
    }

    [[noreturn]] static void throw_at(const SourceLocation& location, const std::string& message) {
        const std::string source = location.str();
        throw std::runtime_error(source.empty() ? message : source + ": " + message);
    }
};

} // namespace fem::io::dsl
