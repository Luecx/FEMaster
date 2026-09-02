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
#include <cmath>
#include <cstdint>
#include <memory>
#include <sstream>
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
            // The first data line contains only the signed term count. Parsing it as
            // signed prevents negative values from wrapping into the unsigned Index type.
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::int64_t>().name("TERMS")
                )
                .bind([ctx](std::int64_t terms) {
                    logging::error(terms >= 2,
                        "EQUATION: at least two terms are required");

                    ctx->remaining    = static_cast<Index>(terms);
                    ctx->first_is_set = false;
                    ctx->equations.clear();
                })
            )
            // Abaqus permits one to four triples per physical data line. Missing tail
            // fields are padded with empty strings so the fixed DSL pattern stays simple.
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<std::string, 12>().name("DATA")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&model, ctx](const std::array<std::string, 12>& data) {
                    logging::error(ctx->remaining > 0,
                        "EQUATION: more term data provided than declared");

                    Index terms_on_line = 0;

                    // Parse each complete target/DOF/coefficient triple and append it
                    // directly to the final equation rows.
                    for (std::size_t i = 0; i < data.size(); i += 3) {
                        if (data[i].empty()) break;

                        logging::error(!data[i + 1].empty() && !data[i + 2].empty(),
                            "EQUATION: incomplete node/DOF/coefficient triple");

                        int       dof         = 0;
                        Precision coefficient = Precision(0);

                        std::istringstream dof_stream        (data[i + 1]);
                        std::istringstream coefficient_stream(data[i + 2]);
                        dof_stream         >> dof;
                        coefficient_stream >> coefficient;

                        logging::error(!dof_stream.fail() && dof_stream.eof() && dof >= 1 && dof <= 6,
                            "EQUATION: DOF must be an integer in [1,6]");
                        logging::error(!coefficient_stream.fail() && coefficient_stream.eof() && std::isfinite(coefficient),
                            "EQUATION: coefficient must be finite");

                        const Dim  equation_dof = static_cast<Dim>(dof - 1);
                        const bool is_set       = model._data->node_sets.has(data[i]);

                        // The first target determines whether this input expands to one
                        // equation or one equation per member of an NSET.
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
                        }
                        // Subsequent NSETs are only legal when the first target was an
                        // NSET, and they must have exactly the same cardinality for pairing.
                        else if (is_set) {
                            logging::error(ctx->first_is_set,
                                "EQUATION: node sets are only valid when the first target is a node set");

                            const auto set = model._data->node_sets.get(data[i]);
                            logging::error(set != nullptr && set->size() == ctx->equations.size(),
                                "EQUATION: node set ", data[i], " has incompatible size");

                            for (std::size_t j = 0; j < ctx->equations.size(); ++j) {
                                ctx->equations[j].entries.push_back({set->at(j), equation_dof, coefficient});
                            }
                        }
                        // A later single node contributes the same term to every row
                        // generated by the first NSET, or to the single scalar equation.
                        else {
                            const ID node = model.compiled_node_id(data[i]);
                            for (auto& equation : ctx->equations) {
                                equation.entries.push_back({node, equation_dof, coefficient});
                            }
                        }

                        ++terms_on_line;
                    }

                    // Validate the physical line before reducing the declared term count.
                    logging::error(terms_on_line > 0,
                        "EQUATION: data line contains no terms");
                    logging::error(terms_on_line <= ctx->remaining,
                        "EQUATION: more terms provided than declared");

                    ctx->remaining -= terms_on_line;
                    if (ctx->remaining != 0) return;

                    // Once all declared terms are present, transfer the completed expanded
                    // rows into the compiled model and clear the temporary parser storage.
                    auto& equations = model._data->equations;
                    for (auto& equation : ctx->equations) {
                        equation.source = constraint::EquationSourceKind::Manual;
                        equations.push_back(std::move(equation));
                    }
                    ctx->equations.clear();
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
