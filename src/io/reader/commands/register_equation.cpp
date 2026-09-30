/**
 * @file register_equation.cpp
 * @brief Registers Abaqus linear constraint equations.
 *
 * Abaqus `EQUATION` records start with the number of terms. Subsequent data
 * lines contain up to four `(node/NSET, dof, coefficient)` triples until the
 * declared number of terms has been collected. Node sets are expanded directly
 * into FEMaster constraint equations while the data lines are consumed.
 *
 * @author Finn Eggers
 * @date 21.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../constraints/types/equation.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../parser_abq.h"

#include <array>
#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader::commands {

/**
 * Registers the Abaqus `EQUATION` grammar in the post-compile analysis pass.
 *
 * The first data line declares the total number of terms. Every following line
 * contributes one to four complete `(target, dof, coefficient)` triples. The
 * first target determines the expansion width: a node creates one constraint
 * equation, while an NSET creates one equation per set entry. Later NSETs are
 * paired by their existing compiled order and later single nodes are reused in
 * every expanded equation.
 *
 * Completed equations are transferred directly into `ModelData::equations`.
 * Parser-local state is retained only until the declared number of terms has
 * been consumed; leaving the command with remaining terms is an input error.
 *
 * @param registry Parser registry receiving the command definition.
 * @param parser Abaqus parser owning the compiled model and analysis state.
 */
void register_equation(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("EQUATION", [&](fem::io::dsl::Command& command) {
        // Temporary state for one Abaqus equation definition. `equations` already
        // contains the final expanded rows; no separate term representation is needed.
        struct Context {
            Index                 remaining    = 0;
            bool                  first_is_set = false;
            constraint::Equations equations;
        };

        auto ctx = std::make_shared<Context>();

        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Define Abaqus linear constraint equations.");

        // EQUATION is interpreted after topology compilation, so all node and NSET
        // references can be resolved immediately into the global assembly namespace.
        command.on_enter([&model, ctx](const fem::io::dsl::Keys&) {
            logging::error(model._data->compiled,
                "EQUATION: requires a compiled model");

            ctx->remaining    = 0;
            ctx->first_is_set = false;
            ctx->equations.clear();
        });

        // A keyword boundary or EOF is only valid after all declared terms were read.
        command.on_exit([ctx](const fem::io::dsl::Keys&) {
            logging::error(ctx->remaining == 0,
                "EQUATION: fewer terms provided than declared");
        });

        command.variant(fem::io::dsl::Variant::make()
            // Every physical line is normalized to twelve string fields. A line with
            // only the first field populated starts a new equation; all other lines
            // contain up to four (target, dof, coefficient) triples.
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<std::string, 12>().name("DATA")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&model, ctx](const std::array<std::string, 12>& data) {
                    if (data[1].empty()) {
                        logging::error(ctx->remaining == 0,
                            "EQUATION: fewer terms provided than declared");

                        const auto terms = std::stoll(data[0]);
                        logging::error(terms >= 2,
                            "EQUATION: at least two terms are required");

                        ctx->remaining    = static_cast<Index>(terms);
                        ctx->first_is_set = false;
                        ctx->equations.clear();
                        return;
                    }

                    logging::error(ctx->remaining > 0,
                        "EQUATION: term data provided before term count");

                    Index terms_on_line = 0;

                    for (std::size_t i = 0; i < data.size() && !data[i].empty(); i += 3) {
                        logging::error(!data[i + 1].empty() && !data[i + 2].empty(),
                            "EQUATION: incomplete node/DOF/coefficient triple");

                        const int       dof         = std::stoi(data[i + 1]);
                        const Precision coefficient = static_cast<Precision>(std::stod(data[i + 2]));

                        logging::error(dof >= 1 && dof <= 6,
                            "EQUATION: DOF must be an integer in [1,6]");

                        const Dim  equation_dof = static_cast<Dim>(dof - 1);
                        const bool is_set       = model._data->node_sets.has(data[i]);

                        if (ctx->equations.empty()) {
                            ctx->first_is_set = is_set;

                            if (is_set) {
                                const auto set = model._data->node_sets.get(data[i]);
                                logging::error(set != nullptr && set->size() > 0,
                                    "EQUATION: first node set is empty");

                                ctx->equations.resize(set->size());
                                for (std::size_t j = 0; j < set->size(); ++j) {
                                    ctx->equations[j].entries.push_back({set->at(j), equation_dof, coefficient});
                                }
                            } else {
                                ctx->equations.resize(1);
                                ctx->equations[0].entries.push_back({
                                    model.compiled_node_id(data[i]), equation_dof, coefficient
                                });
                            }
                        } else if (is_set) {
                            logging::error(ctx->first_is_set,
                                "EQUATION: node sets are only valid when the first target is a node set");

                            const auto set = model._data->node_sets.get(data[i]);
                            logging::error(set != nullptr && set->size() == ctx->equations.size(),
                                "EQUATION: node set ", data[i], " has incompatible size");

                            for (std::size_t j = 0; j < ctx->equations.size(); ++j) {
                                ctx->equations[j].entries.push_back({set->at(j), equation_dof, coefficient});
                            }
                        } else {
                            const ID node = model.compiled_node_id(data[i]);
                            for (auto& equation : ctx->equations) {
                                equation.entries.push_back({node, equation_dof, coefficient});
                            }
                        }

                        ++terms_on_line;
                    }

                    logging::error(terms_on_line <= ctx->remaining,
                        "EQUATION: more terms provided than declared");

                    ctx->remaining -= terms_on_line;
                    if (ctx->remaining != 0) return;

                    for (auto& equation : ctx->equations) {
                        equation.source = constraint::EquationSourceKind::Manual;
                        model._data->equations.push_back(std::move(equation));
                    }
                    ctx->equations.clear();
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
