// register_element.inl — DSL registration for *ELEMENT

#include <array>
#include <string>
#include <utility>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../model/beam/b33.h"
#include "../../../model/model.h"
#include "../../../model/shell/qspt.h"
#include "../../../model/shell/frt_shell_s3.h"
#include "../../../model/shell/frt_shell_s4.h"
#include "../../../model/shell/frt_shell_s6.h"
#include "../../../model/shell/frt_shell_s8.h"
#include "../../../model/shell/s3.h"
#include "../../../model/shell/s4.h"
#include "../../../model/shell/s4_mitc.h"
#include "../../../model/shell/s6.h"
#include "../../../model/shell/s8.h"
#include "../../../model/shell/s8_mitc.h"
#include "../../../model/solid/c3d10.h"
#include "../../../model/solid/c3d15.h"
#include "../../../model/solid/c3d20.h"
#include "../../../model/solid/c3d20r.h"
#include "../../../model/solid/c3d4.h"
#include "../../../model/solid/c3d6.h"
#include "../../../model/solid/c3d8.h"
#include "../../../model/solid/c3d8r.h"
#include "../../../model/truss/truss.h"

namespace fem::io::reader::commands {

template<class Elem, std::size_t N, std::size_t... I>
inline void set_regular_element_impl(model::Model& model,
                                     fem::ID id,
                                     const std::array<fem::ID, N>& nodes,
                                     std::index_sequence<I...>) {
    model.set_element<Elem>(id, nodes[I]...);
}

template<class Elem, std::size_t N>
inline void set_regular_element(model::Model& model, fem::ID id, const std::array<fem::ID, N>& nodes) {
    set_regular_element_impl<Elem, N>(model, id, nodes, std::make_index_sequence<N>{});
}

inline void register_element(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ELEMENT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define finite elements by id and connectivity.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET")
                    .optional("EALL")
                    .doc("Element set to activate while reading")
                .key("TYPE")
                    .required()
                    .doc("Element topology/type")
                    .allowed({
                        "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                        "B33", "T3",
                        "S3", "S4", "MITC4", "S6", "S8", "MITC8", "QSPT",
                        "MITC3FRT", "MITC4FRT", "MITC6FRT", "MITC8FRT"
                    })
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->elem_sets.activate(keys.raw("ELSET"));
        });

#define FEM_ADD_ELEMENT_VARIANT(KEY, ELEM, COUNT, DESC) \
        command.variant(fem::io::dsl::Variant::make() \
            .when(fem::io::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::io::dsl::Segment::make() \
                .range(fem::io::dsl::LineRange{}.min(1)) \
                .pattern(fem::io::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<fem::ID>().name("ID").desc("Element id") \
                    .fixed<fem::ID, COUNT>().name("N").desc(DESC) \
                ) \
                .bind([&model](fem::ID id, const std::array<fem::ID, COUNT>& nodes) { \
                    set_regular_element<model::ELEM>(model, id, nodes); \
                }) \
            ) \
        );

        FEM_ADD_ELEMENT_VARIANT("C3D4", C3D4, 4, "C3D4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D6", C3D6, 6, "C3D6 connectivity (6 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D8", C3D8, 8, "C3D8 connectivity (8 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D8R", C3D8R, 8, "C3D8R connectivity (8 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D10", C3D10, 10, "C3D10 connectivity (10 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D15", C3D15, 15, "C3D15 connectivity (15 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D20", C3D20, 20, "C3D20 connectivity (20 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D20R", C3D20R, 20, "C3D20R connectivity (20 nodes)");
        FEM_ADD_ELEMENT_VARIANT("T3", T3, 2, "T3 connectivity (2 nodes)");
        FEM_ADD_ELEMENT_VARIANT("S3", S3, 3, "S3 connectivity (3 nodes)");
        FEM_ADD_ELEMENT_VARIANT("S4", S4, 4, "S4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("QSPT", QSPT, 4, "QSPT connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("MITC4", MITC4, 4, "MITC4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("MITC3FRT", FRTShellS3, 3, "MITC3FRT connectivity (3 nodes)");
        FEM_ADD_ELEMENT_VARIANT("MITC4FRT", FRTShellS4, 4, "MITC4FRT connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("MITC6FRT", FRTShellS6, 6, "MITC6FRT connectivity (6 nodes)");
        FEM_ADD_ELEMENT_VARIANT("MITC8FRT", FRTShellS8, 8, "MITC8FRT connectivity (8 nodes)");
        FEM_ADD_ELEMENT_VARIANT("S6", S6, 6, "S6 connectivity (6 nodes)");
        FEM_ADD_ELEMENT_VARIANT("S8", S8, 8, "S8 connectivity (8 nodes)");
        FEM_ADD_ELEMENT_VARIANT("MITC8", MITC8, 8, "MITC8 connectivity (8 nodes)");

#undef FEM_ADD_ELEMENT_VARIANT

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"B33"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::ID>().name("ID").desc("Element id")
                    .fixed<fem::ID, 3>()
                        .name("N")
                        .desc("B33 connectivity (2 nodes plus optional orientation node)")
                        .on_missing(fem::ID{-1})
                )
                .bind([&model](fem::ID id, const std::array<fem::ID, 3>& nodes) {
                    const fem::ID orientation = nodes[2];
                    model.set_beam_element<model::B33>(id, orientation, nodes[0], nodes[1]);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"C3D5"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::ID>().name("ID").desc("Element id")
                    .fixed<fem::ID, 5>().name("N").desc("C3D5 connectivity (5 nodes)")
                )
                .bind([&model](fem::ID id, const std::array<fem::ID, 5>& nodes) {
                    model.set_element<model::C3D8>(id,
                                                   nodes[0], nodes[1], nodes[2], nodes[3],
                                                   nodes[4], nodes[4], nodes[4], nodes[4]);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
