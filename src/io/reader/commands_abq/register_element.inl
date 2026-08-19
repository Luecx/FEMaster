/**
 * @file register_element.inl
 * @brief Registers supported Abaqus finite elements.
 *
 * Abaqus element labels are mapped to the corresponding FEMaster beam, truss,
 * shell and solid implementations while connectivity is written into the active
 * semantic Part. Dedicated mappings handle topology differences where an Abaqus
 * label requires adaptation to a supported FEMaster element.
 *
 * Sparse part-local identifiers are retained during parsing. Instance expansion,
 * coordinate placement and dense assembly enumeration remain deferred to
 * `Model::compile()`.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <utility>

#include "../../../model/beam/b33.h"
#include "../../../model/model.h"
#include "../../../model/shell/frt_shell_s3.h"
#include "../../../model/shell/frt_shell_s4.h"
#include "../../../model/shell/frt_shell_s6.h"
#include "../../../model/shell/frt_shell_s8.h"
#include "../../../model/solid/c3d10.h"
#include "../../../model/solid/c3d15.h"
#include "../../../model/solid/c3d20.h"
#include "../../../model/solid/c3d20r.h"
#include "../../../model/solid/c3d4.h"
#include "../../../model/solid/c3d6.h"
#include "../../../model/solid/c3d8.h"
#include "../../../model/solid/c3d8r.h"
#include "../../../model/truss/truss.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

template<class Elem, std::size_t N, std::size_t... I>
inline void set_abq_element_impl(model::Model& model,
                                 ID id,
                                 const std::array<ID, N>& nodes,
                                 std::index_sequence<I...>) {
    model.set_element<Elem>(id, nodes[I]...);
}

template<class Elem, std::size_t N>
inline void set_abq_element(model::Model& model, ID id, const std::array<ID, N>& nodes) {
    set_abq_element_impl<Elem, N>(model, id, nodes, std::make_index_sequence<N>{});
}

inline void register_element(fem::io::dsl::Registry& registry,
                             model::Model& model,
                             std::shared_ptr<bool> assembly_scope) {
    registry.command("ELEMENT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").optional("EALL")
                .key("TYPE").required().allowed({
                    "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                    "B33", "T3D2", "S3", "S3R", "S4", "S4R", "S6", "S6R", "S8", "S8R"
                })
        );
        command.on_enter([&model, assembly_scope](const fem::io::dsl::Keys& keys) {
            logging::error(model._data != nullptr && !model._data->compiled,
                "ELEMENT: elements cannot be added after compile()");
            if (*assembly_scope) model._data->parts.activate(model::Model::DEFAULT_PART_NAME);
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "ELEMENT: no active part is available");
            part->elem_sets.activate(keys.raw("ELSET"));
        });

#define FEM_ABQ_ELEMENT(KEY, ELEM, COUNT) \
        command.variant(fem::io::dsl::Variant::make() \
            .when(fem::io::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::io::dsl::Segment::make() \
                .range(fem::io::dsl::LineRange{}.min(1)) \
                .pattern(fem::io::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<ID>().name("ID") \
                    .fixed<ID, COUNT>().name("N") \
                ) \
                .bind([&model](ID id, const std::array<ID, COUNT>& nodes) { \
                    set_abq_element<model::ELEM>(model, id, nodes); \
                }) \
            ) \
        );

        FEM_ABQ_ELEMENT("C3D4", C3D4, 4);
        FEM_ABQ_ELEMENT("C3D6", C3D6, 6);
        FEM_ABQ_ELEMENT("C3D8", C3D8, 8);
        FEM_ABQ_ELEMENT("C3D8R", C3D8R, 8);
        FEM_ABQ_ELEMENT("C3D10", C3D10, 10);
        FEM_ABQ_ELEMENT("C3D15", C3D15, 15);
        FEM_ABQ_ELEMENT("C3D20", C3D20, 20);
        FEM_ABQ_ELEMENT("C3D20R", C3D20R, 20);
        FEM_ABQ_ELEMENT("B33", B33, 2);
        FEM_ABQ_ELEMENT("T3D2", T3, 2);
        FEM_ABQ_ELEMENT("S3", FRTShellS3, 3);
        FEM_ABQ_ELEMENT("S3R", FRTShellS3, 3);
        FEM_ABQ_ELEMENT("S4", FRTShellS4, 4);
        FEM_ABQ_ELEMENT("S4R", FRTShellS4, 4);
        FEM_ABQ_ELEMENT("S6", FRTShellS6, 6);
        FEM_ABQ_ELEMENT("S6R", FRTShellS6, 6);
        FEM_ABQ_ELEMENT("S8", FRTShellS8, 8);
        FEM_ABQ_ELEMENT("S8R", FRTShellS8, 8);

#undef FEM_ABQ_ELEMENT

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"C3D5"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .one<ID>().name("ID")
                    .fixed<ID, 5>().name("N")
                )
                .bind([&model](ID id, const std::array<ID, 5>& nodes) {
                    model.set_element<model::C3D8>(id,
                        nodes[0], nodes[1], nodes[2], nodes[3],
                        nodes[4], nodes[4], nodes[4], nodes[4]);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
