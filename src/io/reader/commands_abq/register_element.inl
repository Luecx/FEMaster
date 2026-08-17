// register_element.inl — Abaqus *ELEMENT registration
#pragma once

#include <array>
#include <functional>
#include <string>
#include <utility>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../model/beam/b33.h"
#include "../../../model/model.h"
#include "../../../model/shell/s3.h"
#include "../../../model/shell/s4.h"
#include "../../../model/solid/c3d10.h"
#include "../../../model/solid/c3d15.h"
#include "../../../model/solid/c3d20.h"
#include "../../../model/solid/c3d20r.h"
#include "../../../model/solid/c3d4.h"
#include "../../../model/solid/c3d6.h"
#include "../../../model/solid/c3d8.h"
#include "../../../model/solid/c3d8r.h"
#include "../../../model/truss/truss.h"

namespace fem::io::reader::commands_abq {

using ElementCountSink = std::function<void(fem::ID)>;

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

inline void register_element_count(fem::io::dsl::Registry& registry, ElementCountSink sink) {
    registry.command("ELEMENT", [sink = std::move(sink)](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").optional("EALL")
                .key("TYPE").required().allowed({
                    "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                    "B33", "T3D2", "S3", "S4"
                })
        );

#define FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT(KEY, COUNT) \
        command.variant(fem::io::dsl::Variant::make() \
            .when(fem::io::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::io::dsl::Segment::make() \
                .range(fem::io::dsl::LineRange{}.min(1)) \
                .pattern(fem::io::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<fem::ID>().name("ID") \
                    .fixed<fem::ID, COUNT>().name("N") \
                ) \
                .bind([sink](fem::ID id, const std::array<fem::ID, COUNT>&) { sink(id); }) \
            ) \
        );

        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D4", 4);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D5", 5);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D6", 6);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D8", 8);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D8R", 8);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D10", 10);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D15", 15);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D20", 20);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D20R", 20);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("B33", 2);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("T3D2", 2);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S3", 3);
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S4", 4);

#undef FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT
    });
}

inline void register_element(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ELEMENT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define Abaqus finite elements by id and connectivity.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").optional("EALL")
                .key("TYPE").required().allowed({
                    "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                    "B33", "T3D2", "S3", "S4"
                })
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->elem_sets.activate(keys.raw("ELSET"));
        });

#define FEM_ADD_ABQ_ELEMENT_VARIANT(KEY, ELEM, COUNT) \
        command.variant(fem::io::dsl::Variant::make() \
            .when(fem::io::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::io::dsl::Segment::make() \
                .range(fem::io::dsl::LineRange{}.min(1)) \
                .pattern(fem::io::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<fem::ID>().name("ID") \
                    .fixed<fem::ID, COUNT>().name("N") \
                ) \
                .bind([&model](fem::ID id, const std::array<fem::ID, COUNT>& nodes) { \
                    set_regular_element<model::ELEM>(model, id, nodes); \
                }) \
            ) \
        );

        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D4", C3D4, 4);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D6", C3D6, 6);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D8", C3D8, 8);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D8R", C3D8R, 8);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D10", C3D10, 10);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D15", C3D15, 15);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D20", C3D20, 20);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D20R", C3D20R, 20);
        FEM_ADD_ABQ_ELEMENT_VARIANT("T3D2", T3, 2);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S3", S3, 3);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S4", S4, 4);

#undef FEM_ADD_ABQ_ELEMENT_VARIANT

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"C3D5"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::ID>().name("ID")
                    .fixed<fem::ID, 5>().name("N")
                )
                .bind([&model](fem::ID id, const std::array<fem::ID, 5>& nodes) {
                    model.set_element<model::C3D8>(id,
                                                   nodes[0], nodes[1], nodes[2], nodes[3],
                                                   nodes[4], nodes[4], nodes[4], nodes[4]);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"B33"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::ID>().name("ID")
                    .fixed<fem::ID, 2>().name("N")
                )
                .bind([&model](fem::ID id, const std::array<fem::ID, 2>& nodes) {
                    model.set_beam_element<model::B33>(id, fem::ID{-1}, nodes[0], nodes[1]);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
