/**
 * @file deck.h
 * @brief Defines the parsed input-deck tree used by FEMaster readers.
 *
 * `Deck` stores the syntactic result of one complete DSL parse. Each keyword
 * occurrence becomes one `ParsedCommand` with normalized keyword arguments,
 * its parent scope, child commands and the already validated segment invocations
 * belonging to its data block. Physical data lines therefore never become tree
 * nodes; they remain grouped below their owning command as logical invocations.
 *
 * The hierarchy has one real `ROOT` command. Parsed commands own their children
 * directly through `std::unique_ptr`, while non-owning parent pointers provide
 * upward navigation. Node addresses remain stable when child containers grow,
 * so semantic processing can operate directly on command objects without
 * integer node ids or a parallel flat node store.
 *
 * Parsing and semantic execution remain intentionally separated. `DeckParser`
 * determines scopes, variants and record layouts once, while higher-level readers
 * later choose the dependency order through the node-local `children()`,
 * `enter()`, `leave()` and `execute_children()` operations.
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

#include <initializer_list>
#include <memory>
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
 * @brief One keyword occurrence and tree node in a parsed input deck.
 *
 * Each node owns its direct children and stores a non-owning pointer to its
 * parent. The ROOT node is its own parent, so every parsed command always has a
 * valid parent object. Child ownership through `std::unique_ptr` keeps addresses
 * stable while the parser appends further siblings and allows the scope stack to
 * store direct pointers instead of integer ids.
 *
 * The node also owns all logical data-record invocations selected for this
 * keyword. `enter()` executes the command's entry hook followed by those records,
 * while `leave()` executes the exit hook. Children are never traversed implicitly;
 * the semantic reader decides their order explicitly.
 */
class ParsedCommand {
public:
    ParsedCommand() = default;

    ParsedCommand(const ParsedCommand&)            = delete;
    ParsedCommand& operator=(const ParsedCommand&) = delete;
    ParsedCommand(ParsedCommand&&)                 = delete;
    ParsedCommand& operator=(ParsedCommand&&)      = delete;

    // Command identity and tree navigation
    const Command& command() const {
        return *command_;
    }

    const Keys& keys() const {
        return keys_;
    }

    const ParsedCommand& parent() const {
        return *parent_;
    }

    std::vector<const ParsedCommand*> children(const std::string& command = {}) const {
        std::vector<const ParsedCommand*> result;
        for (const auto& child : children_) {
            if (command.empty() || child->command_->name_ == command) {
                result.push_back(child.get());
            }
        }
        return result;
    }

    // Execute this command without implicitly traversing its children
    void enter() const {
        logging::info("Processing command: *", command_->name_);

        if (command_->on_enter_) {
            try {
                command_->on_enter_(parent_info(), keys_);
            } catch (const std::exception& e) {
                throw_at(location_, e.what());
            }
        }

        for (const ParsedInvocation& invocation : invocations_) {
            if (!invocation.segment || !invocation.segment->_invoke) continue;

            try {
                invocation.segment->_invoke(invocation.tokens);
            } catch (const std::exception& e) {
                throw_at(invocation.location, e.what());
            }
        }
    }

    void leave() const {
        if (!command_->on_exit_) return;

        try {
            command_->on_exit_(keys_);
        } catch (const std::exception& e) {
            throw_at(location_, e.what());
        }
    }

    void execute() const {
        enter();
        leave();
    }

    // Execute selected direct children in the requested semantic order
    void execute_children(const std::string& command) const {
        for (const auto& child : children_) {
            if (child->command_->name_ == command) child->execute();
        }
    }

    void execute_children(std::initializer_list<const char*> commands) const {
        for (const char* command : commands) execute_children(command);
    }

    void execute_children() const {
        for (const auto& child : children_) child->execute();
    }

private:
    const Command* command_ = nullptr;
    Keys           keys_;

    ParsedCommand* parent_ = nullptr;

    std::vector<std::unique_ptr<ParsedCommand>> children_;
    std::vector<ParsedInvocation>               invocations_;

    SourceLocation location_;

    friend class Deck;
    friend class DeckParser;

    ParentInfo parent_info() const {
        return ParentInfo{parent_->command_->name_, parent_->keys_};
    }

    [[noreturn]] static void throw_at(const SourceLocation& location, const std::string& message) {
        const std::string source = location.str();
        throw std::runtime_error(source.empty() ? message : source + ": " + message);
    }
};

/**
 * @brief Owns the root of one immutable parsed command tree.
 *
 * The ROOT command is allocated like every other node so its address remains
 * stable when a `Deck` is moved. All remaining nodes are owned recursively by
 * their parent commands. `Deck` therefore only exposes the root; navigation and
 * semantic execution belong directly to `ParsedCommand`.
 */
class Deck {
public:
    Deck()
        : root_(std::make_unique<ParsedCommand>()) {
        root_->command_ = &root_command();
        root_->parent_  = root_.get();
    }

    const ParsedCommand& root() const {
        return *root_;
    }

private:
    std::unique_ptr<ParsedCommand> root_;

    friend class DeckParser;

    static const Command& root_command() {
        static const Command command{"ROOT"};
        return command;
    }
};

} // namespace fem::io::dsl
