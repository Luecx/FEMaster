/**
 * @file register_element.cpp
 * @brief Registers part-local and unqualified assembly finite elements.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <string>
#include <utility>

#include "../../../model/beam/b31.h"
#include "../../../model/beam/b33.h"
#include "../../../model/model.h"
#include "../../../model/shell/frt_shell_s3.h"
#include "../../../model/shell/frt_shell_s4.h"
#include "../../../model/shell/frt_shell_s6.h"
#include "../../../model/shell/frt_shell_s8.h"
#include "../../../model/shell/qspt.h"
#include "../../../model/shell/s3.h"
#include "../../../model/shell/s4.h"
#include "../../../model/shell/s4_mitc.h"
#include "../../../model/shell/s6.h"
#include "../../../model/shell/s8.h"
#include "../../../model/shell/s8_mitc.h"
#include "../../../model/shell/shell_thermal.h"
#include "../../../model/solid/c3d10.h"
#include "../../../model/solid/c3d15.h"
#include "../../../model/solid/c3d20.h"
#include "../../../model/solid/c3d20r.h"
#include "../../../model/solid/c3d4.h"
#include "../../../model/solid/c3d6.h"
#include "../../../model/solid/c3d8.h"
#include "../../../model/solid/c3d8r.h"
#include "../../../model/element/point.h"
#include "../../../model/truss/t3d2_thermal.h"
#include "../../../model/truss/truss.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

template<class Elem, std::size_t N, std::size_t... I>
inline void set_regular_element_impl(model::Model& model,
                                     ID id,
                                     const std::array<ID, N>& nodes,
                                     std::index_sequence<I...>) {
    model.set_element<Elem>(id, nodes[I]...);
}

template<class Elem, std::size_t N>
inline void set_regular_element(model::Model& model, ID id, const std::array<ID, N>& nodes) {
    set_regular_element_impl<Elem, N>(model, id, nodes, std::make_index_sequence<N>{});
}

void register_element(dsl::Registry& registry, model::Model& model) {
    registry.command("ELEMENT", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        command.keyword(
            dsl::KeywordSpec::make()
                .key("ELSET").optional("EALL")
                .key("TYPE").required().allowed({
                    "C3D4", "C3D5", "C3D6", "C3D8", "C3D8R", "C3D10", "C3D15", "C3D20", "C3D20R",
                    "B31", "B33", "T3", "T3D2", "S3", "S4", "MITC4", "S6", "S8", "MITC8", "QSPT",
                    "MITC3FRT", "MITC4FRT", "MITC6FRT", "MITC8FRT", "MASS", "ROTARYI", "SPRING1"
                })
        );

        command.on_enter([&model](const dsl::Keys& keys) {
            logging::error(model._data != nullptr && !model._data->compiled,
                "ELEMENT: elements cannot be added after compile()");

            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "ELEMENT: no active part is available");

            part->elem_sets.activate(keys.raw("ELSET"));
        });

#define FEM_ADD_ELEMENT_VARIANT(KEY, ELEM, COUNT) \
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
                    set_regular_element<model::ELEM>(model, id, nodes); \
                }) \
            ) \
        );

        // Solid formulations already implement ThermalElement on the common solid base.
        FEM_ADD_ELEMENT_VARIANT("C3D4", C3D4, 4);
        FEM_ADD_ELEMENT_VARIANT("C3D6", C3D6, 6);
        FEM_ADD_ELEMENT_VARIANT("C3D8", C3D8, 8);
        FEM_ADD_ELEMENT_VARIANT("C3D8R", C3D8R, 8);
        FEM_ADD_ELEMENT_VARIANT("C3D10", C3D10, 10);
        FEM_ADD_ELEMENT_VARIANT("C3D15", C3D15, 15);
        FEM_ADD_ELEMENT_VARIANT("C3D20", C3D20, 20);
        FEM_ADD_ELEMENT_VARIANT("C3D20R", C3D20R, 20);

        // B31, T3D2, S3 and S4 add thermal conduction while reusing the existing
        // mechanical formulations. The legacy FEMaster-only labels remain unchanged.
        FEM_ADD_ELEMENT_VARIANT("B31", B31, 2);
        FEM_ADD_ELEMENT_VARIANT("B33", B33, 2);
        FEM_ADD_ELEMENT_VARIANT("T3", T3, 2);
        FEM_ADD_ELEMENT_VARIANT("T3D2", ThermalT3D2, 2);

        FEM_ADD_ELEMENT_VARIANT("S3", ThermalS3, 3);
        FEM_ADD_ELEMENT_VARIANT("S4", ThermalS4, 4);
        FEM_ADD_ELEMENT_VARIANT("S6", S6, 6);
        FEM_ADD_ELEMENT_VARIANT("S8", S8, 8);
        FEM_ADD_ELEMENT_VARIANT("QSPT", QSPT, 4);
        FEM_ADD_ELEMENT_VARIANT("MITC4", MITC4, 4);
        FEM_ADD_ELEMENT_VARIANT("MITC3FRT", FRTShellS3, 3);
        FEM_ADD_ELEMENT_VARIANT("MITC4FRT", FRTShellS4, 4);
        FEM_ADD_ELEMENT_VARIANT("MITC6FRT", FRTShellS6, 6);
        FEM_ADD_ELEMENT_VARIANT("MITC8FRT", FRTShellS8, 8);
        FEM_ADD_ELEMENT_VARIANT("MITC8", MITC8, 8);

        FEM_ADD_ELEMENT_VARIANT("MASS", PointElement, 1);
        FEM_ADD_ELEMENT_VARIANT("ROTARYI", PointElement, 1);
        FEM_ADD_ELEMENT_VARIANT("SPRING1", PointElement, 1);

#undef FEM_ADD_ELEMENT_VARIANT

        // C3D5 uses the established degenerate eight-node pyramid mapping and is
        // therefore thermally capable through the C3D8 implementation as well.
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

} // namespace fem::io::reader::commands
