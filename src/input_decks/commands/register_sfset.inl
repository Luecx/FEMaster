#pragma once
/**
 * @file register_sfset.inl
 * @brief Register *SFSET command.
 *
 * Create/activate a surface set and add surface IDs to it.
 *
 * Header:
 *   SFSET|NAME=<set> (default "SFALL")
 *   GENERATE         (flag) turns rows into ranges: start, end[, inc]
 *
 * Data:
 *  - With GENERATE:   start, end[, inc]          (inc defaults to 1)
 *  - Without:         one or more surface IDs per line
 */

#include <array>
#include <string>
#include <vector>
#include <limits>

#include "../../core/logging.h"
#include "../../core/types_num.h"   // fem::ID
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/segment.h"
#include "../../dsl/pattern.h"
#include "../../dsl/variant.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_sfset(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("SFSET", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Create or activate a surface set and add surface IDs (lists or ranges).");

        // --- Header keywords ---
        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("SFSET").alternative("NAME").optional("SFALL")
                    .doc("Surface set name (default: SFALL).")
                .flag("GENERATE")
                    .doc("Interpret data lines as ranges start,end[,inc].")
        );

        // --- On enter: activate set ---
        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string set_name = keys.get<std::string>("SFSET");
            model._data->surface_sets.activate(set_name);
        });

        // Sentinel used for “explicit IDs” variant to pad shorter lines.
        const fem::ID missing_id = std::numeric_limits<fem::ID>::min();

        // --- Variant A: GENERATE with 3 values: start,end,inc ---
        // Boolean-aware admission: matches when GENERATE is present and truthy
        // (flag without value OR 1/ON/YES/TRUE per Keys::get<bool> semantics).
        command.variant(
            fem::dsl::Variant::make()
                .rank(20)
                .when(fem::dsl::Condition::key_true("GENERATE"))
                .doc("Range rows with explicit increment: start, end, inc.")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(0))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::ID, 3>().name("RANGE3").desc("start, end, inc")
                        )
                        .bind([&model](const std::array<fem::ID, 3>& r) {
                            fem::ID start = r[0], end = r[1], inc = r[2];
                            if (inc == 0) {
                                logging::error(false, "SFSET/GENERATE: increment must not be zero.");
                                return;
                            }
                            if ((inc > 0 && start > end) || (inc < 0 && start < end)) {
                                logging::warning(false, "SFSET/GENERATE: range likely empty with given increment.");
                            }
                            for (fem::ID sid = start; (inc > 0) ? (sid <= end) : (sid >= end); sid += inc) {
                                model._data->surface_sets.get()->add(sid);
                            }
                        })
                )
        );

        // --- Variant B: GENERATE with 2 values: start,end (inc=1 default) ---
        // Also gated by GENERATE being truthy; lower rank so A is tried first.
        command.variant(
            fem::dsl::Variant::make()
                .rank(10)
                .when(fem::dsl::Condition::key_true("GENERATE"))
                .doc("Range rows with implicit increment: start, end (inc=1).")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(0))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::ID, 2>().name("RANGE2").desc("start, end (inc=1)")
                        )
                        .bind([&model](const std::array<fem::ID, 2>& r) {
                            const fem::ID start = r[0], end = r[1];
                            const fem::ID inc   = (start <= end) ? fem::ID{1} : fem::ID{-1};
                            for (fem::ID sid = start; (inc > 0) ? (sid <= end) : (sid >= end); sid += inc) {
                                model._data->surface_sets.get()->add(sid);
                            }
                        })
                )
        );

        // --- Variant C: explicit surface IDs (no GENERATE), multiple per line ---
        // Matches when GENERATE is missing OR explicitly false (0/OFF/NO/FALSE).
        command.variant(
            fem::dsl::Variant::make()
                .when(fem::dsl::Condition::any_of({
                    fem::dsl::Condition::negate(fem::dsl::Condition::key_present("GENERATE")),
                    fem::dsl::Condition::key_false("GENERATE")
                }))
                .doc("Explicit surface IDs (one or more per line).")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(0))
                        .pattern(
                            // Use a bounded-width fixed block with padding so we can accept 0..32 IDs per line.
                            fem::dsl::Pattern::make()
                                .fixed<fem::ID, 32>()
                                    .name("IDS")
                                    .desc("Surface IDs (up to 32 per line).")
                                    .on_missing(missing_id)
                                    .on_empty(missing_id)
                        )
                        .bind([&model, missing_id](const std::array<fem::ID, 32>& ids) {
                            auto set = model._data->surface_sets.get();
                            if (!set) return;
                            for (const fem::ID sid : ids) {
                                if (sid == missing_id) continue;
                                set->add(sid);
                            }
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
