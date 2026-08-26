/**
 * @file engine.h
 * @brief Executes one parsed DSL deck in its original hierarchical order.
 *
 * `Engine` is the compact convenience path for callers that do not need a
 * custom semantic dependency order. It parses the source once through
 * `DeckParser` and then recursively executes every stored command occurrence in
 * source order, preserving entry/exit scope semantics.
 *
 * FEMaster readers use `DeckParser` directly because model construction requires
 * an explicit order around `Model::compile()`. The engine intentionally contains
 * no command activation modes or parser stages.
 *
 * @see DeckParser
 * @see Deck
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "deck_parser.h"
#include "file.h"
#include "registry.h"

namespace fem::io::dsl {

/**
 * @brief Convenience parser/executor for source-ordered DSL processing.
 *
 * Commands with children remain entered while those children are executed and
 * are left afterwards. This reproduces ordinary nested-scope behavior while the
 * syntax itself is still parsed and validated only once.
 */
class Engine {
public:
    explicit Engine(const Registry& registry) : registry_(registry) {}

    /**
     * @brief Parses the complete file once and executes all parsed commands.
     */
    void run(File& file) const {
        const Deck deck = DeckParser(registry_).parse(file);
        for (const Deck::NodeId root : deck.roots()) execute_subtree(deck, root);
    }

private:
    const Registry& registry_;

    static void execute_subtree(const Deck& deck, Deck::NodeId node) {
        deck.enter(node);
        for (const Deck::NodeId child : deck.children(node)) execute_subtree(deck, child);
        deck.leave(node);
    }
};

} // namespace fem::io::dsl
