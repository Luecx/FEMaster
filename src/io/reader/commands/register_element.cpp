/**
 * @file register_element.cpp
 * @brief Registers part-local and unqualified assembly finite elements.
 *
 * The `ELEMENT` command dispatches supported FEMaster type names to their
 * concrete beam, truss, shell and solid element classes. Connectivity is stored
 * in the active semantic Part before compilation; unqualified root definitions
 * use the model's default part.
 *
 * One-node point-element labels `MASS`, `ROTARYI` and `SPRING1` share the
 * zero-dimensional `model::PointElement` implementation. Their physical mass,
 * rotary-inertia or ground-spring contribution is supplied later by the
 * corresponding element-property command.
 *
 * The registration preserves sparse user identifiers and natural connectivity.
 * Instance expansion and dense global enumeration are deliberately deferred to
 * `Model::compile()`.
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <initializer_list>
#include <string>
#include <tuple>
#include <type_traits>

#include "../../../model/beam/b31.h"
#include "../../../model/beam/b33.h"
#include "../../../model/model.h"
#include "../../../model/shell/frt_shell_s3.h"
#include "../../../model/shell/frt_shell_s4.h"
#include "../../../model/shell/frt_shell_s6.h"
#include "../../../model/shell/frt_shell_s8.h"
#include "../../../model/shell/qspt.h"
#include "../../../model/solid/c3d10.h"
#include "../../../model/solid/c3d13.h"
#include "../../../model/solid/c3d15.h"
#include "../../../model/solid/c3d20.h"
#include "../../../model/solid/c3d20r.h"
#include "../../../model/solid/c3d4.h"
#include "../../../model/solid/c3d5.h"
#include "../../../model/solid/c3d6.h"
#include "../../../model/solid/c3d8.h"
#include "../../../model/solid/c3d8i.h"
#include "../../../model/solid/c3d8r.h"
#include "../../../model/element/point.h"
#include "../../../model/truss/truss.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers sparse finite-element connectivity before model compilation.
 */
void register_element(dsl::Registry& registry, model::Model& model) {
    registry.command("ELEMENT", [&](dsl::Command& command) {
        // Elements may be defined directly, inside a Part or inside the assembly scope
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        // Define the destination element set and concrete element formulation
        command.keyword(dsl::KeywordSpec::make()
            .key("ELSET").optional("EALL")
            .key("TYPE").required().allowed({
                "C3D4"     , "C3D5"     , "C3D6"     , "C3D8"     ,
                "C3D8I"    , "C3D8R"    , "C3D10"    , "C3D13"    ,
                "C3D15"    , "C3D20"    , "C3D20R"   ,
                "B31"      , "B33"      , "T3"       , "T3D2"     ,
                "S3"       , "S3R"      , "S4"       , "S4R"      ,
                "S6"       , "S6R"      , "S8"       , "S8R"      ,
                "MITC4"    , "MITC8"    , "MITC3FRT" , "MITC4FRT" ,
                "MITC6FRT" , "MITC8FRT" , "QSPT"     ,
                "MASS"     , "ROTARYI"  , "SPRING1"
            })
        );

        command.on_enter([&model](const dsl::Keys& keys) {
            const auto part = model._data->parts.get();

            logging::error(!model._data->compiled,
                "ELEMENT: elements cannot be added after compile()");
            logging::error(part != nullptr,
                "ELEMENT: no active part is available");

            part->elem_sets.activate(keys.raw("ELSET"));
        });

        // Local generic lambda: the array tag encodes type and connectivity size.
        auto register_variant = [&](auto tag, std::initializer_list<std::string> types) {
            using Element = std::remove_pointer_t<typename decltype(tag)::value_type>;
            constexpr std::size_t N = std::tuple_size_v<decltype(tag)>;

            command.variant(dsl::Variant::make()
                .when(dsl::Condition::key_equals("TYPE", types))
                .segment(dsl::Segment::make()
                    .range(dsl::LineRange{}.min(1))
                    .pattern(dsl::Pattern::make()
                        .allow_multiline()
                        .one<ID>().name("ID")
                        .fixed<ID, N>().name("N")
                    )
                    .bind([&model](ID id, const std::array<ID, N>& nodes) {
                        std::apply([&](auto... node_ids) {
                            model.set_element<Element>(id, node_ids...);
                        }, nodes);
                    })
                )
            );
        };

        // Solids
        register_variant(std::array<model::C3D4*,    4>{}, {"C3D4"});
        register_variant(std::array<model::C3D5*,    5>{}, {"C3D5"});
        register_variant(std::array<model::C3D6*,    6>{}, {"C3D6"});
        register_variant(std::array<model::C3D8*,    8>{}, {"C3D8"});
        register_variant(std::array<model::C3D8I*,   8>{}, {"C3D8I"});
        register_variant(std::array<model::C3D8R*,   8>{}, {"C3D8R"});
        register_variant(std::array<model::C3D10*,  10>{}, {"C3D10"});
        register_variant(std::array<model::C3D13*,  13>{}, {"C3D13"});
        register_variant(std::array<model::C3D15*,  15>{}, {"C3D15"});
        register_variant(std::array<model::C3D20*,  20>{}, {"C3D20"});
        register_variant(std::array<model::C3D20R*, 20>{}, {"C3D20R"});

        // Beams and trusses
        register_variant(std::array<model::B31*, 2>{}, {"B31"});
        register_variant(std::array<model::B33*, 2>{}, {"B33"});
        register_variant(std::array<model::T3*,  2>{}, {"T3", "T3D2"});

        // Shells
        register_variant(std::array<model::FRTShellS3*, 3>{}, {"S3", "S3R", "MITC3FRT"});
        register_variant(std::array<model::FRTShellS4*, 4>{}, {"S4", "S4R", "MITC4", "MITC4FRT"});
        register_variant(std::array<model::FRTShellS6*, 6>{}, {"S6", "S6R", "MITC6FRT"});
        register_variant(std::array<model::FRTShellS8*, 8>{}, {"S8", "S8R", "MITC8", "MITC8FRT"});
        register_variant(std::array<model::QSPT*,       4>{}, {"QSPT"});

        // Point elements
        register_variant(std::array<model::PointElement*, 1>{}, {"MASS", "ROTARYI", "SPRING1"});
    });
}

} // namespace fem::io::reader::commands
