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
 * Parsing and semantic execution are intentionally separated. `DeckParser`
 * determines scopes, variants and record layouts once, while the higher-level
 * FEMaster reader later chooses the dependency order in which stored commands
 * are executed. This allows Part-local topology to be built before
 * `Model::compile()` and assembly-level records afterwards without replaying or
 * re-tokenizing the input file.
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
#include <limits>
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
 * physical data lines never become tree nodes. Closing keywords such as
 * `ENDPART` and `ENDSTEP` are retained as ordinary children of the scope they
 * terminate so any command-specific callback remains available during later
 * semantic processing.
 */
struct ParsedCommand {
    using NodeId = std::size_t;

    const Command* command = nullptr;
    Keys           keys;

    NodeId parent = std::numeric_limits<NodeId>::max();

    std::vector<NodeId>           children;
    std::vector<ParsedInvocation> invocations;

    SourceLocation location;
};

/**
 * @brief Immutable parsed deck with explicit semantic execution operations.
 *
 * `Deck` owns the ordered command occurrences produced by `DeckParser`. Queries
 * return stable integer node ids, while `enter()`, `leave()` and `execute()`
 * invoke the callbacks already registered on the corresponding DSL command.
 * The reader can therefore expose its dependency order directly, for example by
 * entering one Part, executing NODE/ELEMENT/NSET children in a chosen order and
 * leaving that Part before moving to the next one.
 */
class Deck {
public:
    using NodeId = ParsedCommand::NodeId;
    static constexpr NodeId ROOT = std::numeric_limits<NodeId>::max();

    // Access parsed command occurrences and filtered root/child ranges
    const ParsedCommand& node(NodeId id) const {
        return nodes_.at(id);
    }

    std::vector<NodeId> roots(const std::string& command = {}) const {
        std::vector<NodeId> result;
        for (const NodeId id : roots_) {
            if (command.empty() || nodes_[id].command->name_ == command) {
                result.push_back(id);
            }
        }
        return result;
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

    // Execute selected root or direct-child keywords in the requested order
    void execute_roots(const std::string& command) const {
        for (const NodeId id : roots(command)) execute(id);
    }

    void execute_roots(std::initializer_list<const char*> commands) const {
        for (const char* command : commands) execute_roots(command);
    }

    void execute_children(NodeId parent, const std::string& command) const {
        for (const NodeId id : children(parent, command)) execute(id);
    }

    void execute_children(NodeId parent, std::initializer_list<const char*> commands) const {
        for (const char* command : commands) execute_children(parent, command);
    }

    void execute_children(NodeId parent) const {
        for (const NodeId id : children(parent)) execute(id);
    }

private:
    std::vector<ParsedCommand> nodes_;
    std::vector<NodeId>        roots_;

    friend class DeckParser;

    ParentInfo parent_info(const ParsedCommand& current) const {
        if (current.parent == ROOT) return ParentInfo{"ROOT", Keys{}};

        const ParsedCommand& parent = nodes_.at(current.parent);
        return ParentInfo{parent.command->name_, parent.keys};
    }

    [[noreturn]] static void throw_at(const SourceLocation& location, const std::string& message) {
        const std::string source = location.str();
        throw std::runtime_error(source.empty() ? message : source + ": " + message);
    }
};

} // namespace fem::io::dsl
