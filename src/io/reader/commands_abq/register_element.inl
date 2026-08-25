/**
 * @file register_element.inl
 * @brief Registers supported Abaqus finite elements.
 *
 * Abaqus element labels are mapped to the corresponding FEMaster beam, truss,
 * shell and solid implementations while connectivity is written into the active
 * semantic Part. Dedicated mappings handle topology differences where an Abaqus
 * label requires adaptation to a supported FEMaster element.
 *
 * Abaqus `MASS` is represented by the one-node zero-dimensional `PointElement`.
 * The point element preserves element identity, ELSET membership and instance
 * expansion but contributes no matrices or active DOFs itself. A later `*MASS`
 * property reads its node connectivity and creates a separate nodal
 * `feature::PointMass`.
 *
 * Sparse part-local identifiers are retained during parsing. Instance expansion,
 * coordinate placement and dense assembly enumeration remain deferred to
 * `Model::compile()`.
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include <array>
#include <utility>

#include "../../../model/beam/b33.h"
#include "../../../model/element/point.h"
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

namespace dsl = fem::io::dsl;

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

/**
 * @brief Registers Abaqus element connectivity before model compilation.
 */
inline void register_element(dsl::Registry& registry, model::Model& model) {
    registry.command("ELEMENT", [&](dsl::Command& command) {
        // Abaqus elements may be defined directly, in a Part or in the assembly scope
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        // Define the destination set and supported Abaqus element labels
        command.keyword(
            dsl::KeywordSpec::make()
                .key("ELSET").optional("EALL")
                .key("TYPE").required().allowed({
                    "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                    "B33", "T3D2", "S3", "S3R", "S4", "S4R", "S6", "S6R", "S8", "S8R", "MASS"
                })
        );

        // Prepare the active element set before processing connectivity rows
        command.on_enter([&model](const dsl::Keys& keys) {
            logging::error(model._data != nullptr && !model._data->compiled,
                "ELEMENT: elements cannot be added after compile()");

            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "ELEMENT: no active part is available");

            part->elem_sets.activate(keys.raw("ELSET"));
        });

#define FEM_ABQ_ELEMENT(KEY, ELEM, COUNT) \
        command.variant(dsl::Variant::make() \
            .when(dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(dsl::Segment::make() \
                .range(dsl::LineRange{}.min(1)) \
                .pattern(dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<ID>().name("ID") \
                    .fixed<ID, COUNT>().name("N") \
                ) \
                .bind([&model](ID id, const std::array<ID, COUNT>& nodes) { \
                    set_abq_element<model::ELEM>(model, id, nodes); \
                }) \
            ) \
        );

        // Register direct Abaqus-to-FEMaster topology mappings
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
        FEM_ABQ_ELEMENT("MASS", PointElement, 1);

#undef FEM_ABQ_ELEMENT

        // Expand Abaqus C3D5 into the supported degenerate C3D8 topology
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::key_equals("TYPE", {"C3D5"}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
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
