// register_coupling.inl — registers *COUPLING

#include <array>
#include <stdexcept>
#include <string>

#include "../../constraints/coupling.h"
#include "../../core/types_eig.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_coupling(fem::dsl::Registry& registry, model::Model& model) {
    std::string master;
    std::string slave;
    std::string surface;
    std::string type;

    registry.command("COUPLING", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Create kinematic/structural couplings between sets.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MASTER").required().doc("Master node set")
                .key("TYPE").required().allowed({"KINEMATIC", "STRUCTURAL"})
                .key("SLAVE").optional().doc("Slave node set")
                .key("SFSET").optional().doc("Slave surface set")
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            master = keys.raw("MASTER");
            type   = keys.raw("TYPE");
            slave  = keys.has("SLAVE") ? keys.raw("SLAVE") : std::string{};
            surface = keys.has("SFSET") ? keys.raw("SFSET") : std::string{};
            if (slave.empty() == surface.empty()) {
                throw std::runtime_error("COUPLING requires exactly one of SLAVE or SFSET");
            }
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .fixed<fem::Precision, 6>().name("DOF").desc("Coupled degrees of freedom (1=on,0=off)")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&](const std::array<fem::Precision, 6>& dofs_raw) {
                    fem::Dofs mask;
                    for (int i = 0; i < 6; ++i) {
                        mask(i) = dofs_raw[i] > fem::Precision{0};
                    }

                    const bool is_surface = !surface.empty();
                    const std::string& slave_ref = is_surface ? surface : slave;

                    constraint::CouplingType ctype;
                    if (type == "KINEMATIC") {
                        ctype = constraint::CouplingType::KINEMATIC;
                    } else if (type == "STRUCTURAL") {
                        ctype = constraint::CouplingType::STRUCTURAL;
                    } else {
                        throw std::runtime_error("Unsupported COUPLING TYPE=" + type);
                    }

                    model.add_coupling(master, slave_ref, mask, ctype, is_surface);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands

