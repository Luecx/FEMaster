/**
 * @file register_element.inl
 * @brief Registers Abaqus *ELEMENT parsing for allocation and topology construction.
 *
 * The Abaqus element command is handled separately from the native FEMaster
 * command because Abaqus element labels do not map one-to-one to the internal
 * FEMaster element names. The count registration reads the same connectivity
 * layouts without constructing elements, while the topology registration maps
 * supported Abaqus labels directly onto the corresponding FEMaster element
 * implementations.
 *
 * Conventional shell labels are intentionally mapped to the finite-rotation
 * FRT shell family. Full- and reduced-integration Abaqus labels of the same
 * topology therefore create the same FEMaster FRT shell implementation.
 *
 * @see model::FRTShellS3
 * @see model::FRTShellS4
 * @see model::FRTShellS6
 * @see model::FRTShellS8
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
#include <functional>
#include <utility>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
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

namespace fem::io::reader::commands_abq {

/**
 * Registers the allocation-pass form of Abaqus `*ELEMENT`.
 *
 * The command validates the supported Abaqus element label and consumes the
 * complete connectivity record, but forwards only the element identifier to
 * the supplied sink. No model topology is modified during this pass.
 *
 * Shell labels with and without the reduced-integration suffix use identical
 * connectivity sizes because both variants are mapped onto the same FRT shell
 * topology during the later topology pass.
 *
 * @param registry Stage-local DSL registry.
 * @param sink Callback receiving every parsed element identifier.
 */
inline void register_element_count(fem::io::dsl::Registry& registry, std::function<void(fem::ID)> sink) {
    registry.command("ELEMENT", [sink = std::move(sink)](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Count supported Abaqus element identifiers without constructing topology.");

        // Accept only the Abaqus element labels currently mapped by FEMaster
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET")
                    .optional("EALL")
                    .doc("Element set populated by the element definitions")
                .key("TYPE")
                    .required()
                    .doc("Abaqus element type")
                    .allowed({
                        "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                        "B33", "T3D2",
                        "S3", "S3R", "S4", "S4R", "S6", "S6R", "S8", "S8R"
                    })
        );

        // Parse the complete connectivity so malformed element records fail in
        // the count pass before model storage is allocated
#define FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT(KEY, COUNT, DESC) \
        command.variant(fem::io::dsl::Variant::make() \
            .when(fem::io::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::io::dsl::Segment::make() \
                .range(fem::io::dsl::LineRange{}.min(1)) \
                .pattern(fem::io::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<fem::ID>().name("ID").desc("Element identifier") \
                    .fixed<fem::ID, COUNT>().name("N").desc(DESC) \
                ) \
                .bind([sink](fem::ID id, const std::array<fem::ID, COUNT>&) { \
                    sink(id); \
                }) \
            ) \
        );

        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D4"  ,  4, "C3D4 connectivity (4 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D5"  ,  5, "C3D5 connectivity (5 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D6"  ,  6, "C3D6 connectivity (6 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D8"  ,  8, "C3D8 connectivity (8 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D8R" ,  8, "C3D8R connectivity (8 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D10" , 10, "C3D10 connectivity (10 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D15" , 15, "C3D15 connectivity (15 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D20" , 20, "C3D20 connectivity (20 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("C3D20R", 20, "C3D20R connectivity (20 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("B33"   ,  2, "B33 connectivity (2 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("T3D2"  ,  2, "T3D2 connectivity (2 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S3"    ,  3, "S3 connectivity (3 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S3R"   ,  3, "S3R connectivity (3 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S4"    ,  4, "S4 connectivity (4 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S4R"   ,  4, "S4R connectivity (4 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S6"    ,  6, "S6 connectivity (6 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S6R"   ,  6, "S6R connectivity (6 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S8"    ,  8, "S8 connectivity (8 nodes)");
        FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT("S8R"   ,  8, "S8R connectivity (8 nodes)");

#undef FEM_ADD_ABQ_ELEMENT_COUNT_VARIANT
    });
}

/**
 * Registers topology construction for supported Abaqus `*ELEMENT` definitions.
 *
 * Element identifiers and connectivity are written directly into the active
 * FEMaster model. Abaqus `T3D2` is represented by the native two-node truss,
 * while the supported conventional shell labels are represented by the FRT
 * shell family. The reduced-integration suffix does not select a different
 * FEMaster shell formulation.
 *
 * `C3D5` retains the existing FEMaster pyramid representation by degenerating a
 * C3D8 and repeating the apex node. Abaqus `B33` uses the two connectivity nodes
 * directly and leaves the optional FEMaster orientation-node identifier unset.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the constructed topology.
 */
inline void register_element(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ELEMENT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define supported Abaqus finite elements by identifier and connectivity.");

        // Validate the supported Abaqus labels and activate the optional ELSET
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET")
                    .optional("EALL")
                    .doc("Element set populated by the element definitions")
                .key("TYPE")
                    .required()
                    .doc("Abaqus element type")
                    .allowed({
                        "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                        "B33", "T3D2",
                        "S3", "S3R", "S4", "S4R", "S6", "S6R", "S8", "S8R"
                    })
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->elem_sets.activate(keys.raw("ELSET"));
        });

        // Construct element types whose Abaqus connectivity maps directly onto
        // one FEMaster element implementation
#define FEM_ADD_ABQ_ELEMENT_VARIANT(KEY, ELEM, COUNT, DESC, ...) \
        command.variant(fem::io::dsl::Variant::make() \
            .when(fem::io::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::io::dsl::Segment::make() \
                .range(fem::io::dsl::LineRange{}.min(1)) \
                .pattern(fem::io::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<fem::ID>().name("ID").desc("Element identifier") \
                    .fixed<fem::ID, COUNT>().name("N").desc(DESC) \
                ) \
                .bind([&model](fem::ID id, const std::array<fem::ID, COUNT>& nodes) { \
                    model.set_element<model::ELEM>(id, __VA_ARGS__); \
                }) \
            ) \
        );

        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D4"  , C3D4      ,  4, "C3D4 connectivity (4 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D6"  , C3D6      ,  6, "C3D6 connectivity (6 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D8"  , C3D8      ,  8, "C3D8 connectivity (8 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D8R" , C3D8R     ,  8, "C3D8R connectivity (8 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D10" , C3D10     , 10, "C3D10 connectivity (10 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7], nodes[8], nodes[9]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D15" , C3D15     , 15, "C3D15 connectivity (15 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7],
                                    nodes[8], nodes[9], nodes[10], nodes[11], nodes[12], nodes[13], nodes[14]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D20" , C3D20     , 20, "C3D20 connectivity (20 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7],
                                    nodes[8], nodes[9], nodes[10], nodes[11], nodes[12], nodes[13], nodes[14], nodes[15],
                                    nodes[16], nodes[17], nodes[18], nodes[19]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("C3D20R", C3D20R    , 20, "C3D20R connectivity (20 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7],
                                    nodes[8], nodes[9], nodes[10], nodes[11], nodes[12], nodes[13], nodes[14], nodes[15],
                                    nodes[16], nodes[17], nodes[18], nodes[19]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("T3D2"  , T3        ,  2, "T3D2 connectivity (2 nodes)", nodes[0], nodes[1]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S3"    , FRTShellS3,  3, "S3 connectivity (3 nodes)", nodes[0], nodes[1], nodes[2]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S3R"   , FRTShellS3,  3, "S3R connectivity (3 nodes)", nodes[0], nodes[1], nodes[2]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S4"    , FRTShellS4,  4, "S4 connectivity (4 nodes)", nodes[0], nodes[1], nodes[2], nodes[3]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S4R"   , FRTShellS4,  4, "S4R connectivity (4 nodes)", nodes[0], nodes[1], nodes[2], nodes[3]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S6"    , FRTShellS6,  6, "S6 connectivity (6 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S6R"   , FRTShellS6,  6, "S6R connectivity (6 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S8"    , FRTShellS8,  8, "S8 connectivity (8 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
        FEM_ADD_ABQ_ELEMENT_VARIANT("S8R"   , FRTShellS8,  8, "S8R connectivity (8 nodes)",
                                    nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);

#undef FEM_ADD_ABQ_ELEMENT_VARIANT

        // Represent the five-node pyramid with the established degenerated C3D8
        // connectivity used by the native FEMaster reader
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"C3D5"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::ID>().name("ID").desc("Element identifier")
                    .fixed<fem::ID, 5>().name("N").desc("C3D5 connectivity (5 nodes)")
                )
                .bind([&model](fem::ID id, const std::array<fem::ID, 5>& nodes) {
                    model.set_element<model::C3D8>(
                        id,
                        nodes[0], nodes[1], nodes[2], nodes[3],
                        nodes[4], nodes[4], nodes[4], nodes[4]
                    );
                })
            )
        );

        // Abaqus B33 has no explicit FEMaster orientation-node entry in the
        // connectivity record, so leave the optional identifier unset
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"B33"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::ID>().name("ID").desc("Element identifier")
                    .fixed<fem::ID, 2>().name("N").desc("B33 connectivity (2 nodes)")
                )
                .bind([&model](fem::ID id, const std::array<fem::ID, 2>& nodes) {
                    model.set_beam_element<model::B33>(id, fem::ID{-1}, nodes[0], nodes[1]);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
