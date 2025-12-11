#pragma once
/**
 * @file register_surface.inl
 * @brief Register *SURFACE command.
 *
 * Define surfaces from element faces or element sets.
 * Header:
 *   SFSET|NAME=<set> (optional, default "SFALL")
 *
 * Data (repeatable lines):
 *   1) ID, ELEM_ID, SIDE
 *      - SIDE accepts "S#" or integer
 *      - calls model.set_surface(id, elem_id, side)
 *   2) TARGET, SIDE
 *      - TARGET = ELSET name or single element id (int)
 *      - if ELSET exists: model.set_surface(set_name, side)
 *      - else if integer: model.set_surface(-1, elem_id, side)
 *      - else: error
 */

#include <stdexcept>
#include <string>

#include "../../core/logging.h"
#include "../../core/types_num.h"    // fem::ID
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/segment.h"
#include "../../dsl/pattern.h"
#include "../../dsl/variant.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_surface(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("SURFACE", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define surfaces from element faces (by element id) or entire element sets.\n"
            "TYPE must be ELEMENT. Use SFSET/NAME to select the active surface set."
        );

        // --- Header keywords ---
        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("SFSET").alternative("NAME").optional("SFALL")
                    .doc("Set name to activate/create (default: SFALL).\n"
                         "Elements with 2D faces contribute surfaces; beam elements contribute 1D lines under the same set name.")
        );

        // --- On enter: validate + activate sfset ---
        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& sfset = keys.get<std::string>("SFSET"); // resolves NAME alias
            // Activate both surface and line sets under the same name so beam elements
            // can contribute their 1D geometry to a parallel line-set with identical key.
            auto s = model._data->surface_sets.activate(sfset);
            auto l = model._data->line_sets.activate(sfset);
            if (s) { s->sorted(true).duplicates(false); }
            if (l) { l->sorted(true).duplicates(false); }
        });

        // Helper: parse SIDE tokens like "S1", integer, and SPOS/SNEG (case-insensitive)
        auto parse_side = [](const std::string& side_tok) -> int {
            if (side_tok.empty()) return 1;
            std::string s = side_tok;
            for (char& c : s) c = static_cast<char>(std::toupper(c));
            if (s == "SPOS") return 1; // positive normal (shells)
            if (s == "SNEG") return 2; // negative normal (shells use !=1 to flip)
            if (s[0] == 'S' && s.size() > 1) {
                return std::stoi(s.substr(1));
            }
            return std::stoi(s);
        };

        // --- Variant 1: ID, ELEM_ID, SIDE ---
        command.variant(
            fem::dsl::Variant::make()
                .doc("Three fields per line: ID, ELEM_ID, SIDE (SIDE accepts 'S#' or integer).")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .one<fem::ID>().name("ID").desc("Surface id")
                                .one<fem::ID>().name("ELEM_ID").desc("Element id")
                                .one<std::string>().name("SIDE").desc("Face id, e.g. S1 or 1")
                        )
                        .bind([&model, parse_side](fem::ID id, fem::ID elem_id, const std::string& side_tok) {
                            int side = parse_side(side_tok);
                            model.set_surface(id, elem_id, side);
                        })
                )
        );

        // --- Variant 2: TARGET, SIDE (id := -1 if TARGET is elem id) ---
        command.variant(
            fem::dsl::Variant::make()
                .doc("Two fields per line: TARGET, SIDE (TARGET = ELSET name or single element id).")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .one<std::string>().name("TARGET").desc("ELSET name or element id (int)")
                                .one<std::string>().name("SIDE").desc("Face id, e.g. S1 or 1")
                        )
                        .bind([&model, parse_side](const std::string& target, const std::string& side_tok) {
                            int side = parse_side(side_tok);

                            // First: ELSET name?
                            if (model._data->elem_sets.has(target)) {
                                model.set_surface(target, side);
                                return;
                            }

                            // Else: try single element id
                            try {
                                const fem::ID elem_id = static_cast<fem::ID>(std::stoi(target));
                                model.set_surface(-1, elem_id, side); // id := -1 by spec
                            } catch (const std::exception&) {
                                throw std::runtime_error(
                                    "SURFACE target '" + target +
                                    "' is neither an existing element set nor a valid element id."
                                );
                            }
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
