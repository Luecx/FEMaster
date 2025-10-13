// register_element.inl — DSL registration for *ELEMENT

#include <array>
#include <stdexcept>
#include <string>
#include <tuple>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include "../../model/beam/b33.h"
#include "../../model/model.h"
#include "../../model/pointelem/point.h"
#include "../../model/shell/s3.h"
#include "../../model/shell/s4.h"
#include "../../model/shell/s4_mitc.h"
#include "../../model/shell/s6.h"
#include "../../model/shell/s8.h"
#include "../../model/solid/c3d10.h"
#include "../../model/solid/c3d15.h"
#include "../../model/solid/c3d20.h"
#include "../../model/solid/c3d20r.h"
#include "../../model/solid/c3d4.h"
#include "../../model/solid/c3d6.h"
#include "../../model/solid/c3d8.h"
#include "../../model/truss/truss.h"

namespace fem::input_decks::commands {

namespace detail {

template<class Elem, std::size_t NodeCount>
void set_element_from_array(fem::model::Model& model, fem::ID id, const std::array<fem::ID, NodeCount>& nodes) {
    std::apply([&](auto... ns) {
        model.set_element<Elem>(id, ns...);
    }, nodes);
}

} // namespace detail

inline void register_element(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("ELEMENT", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Define finite elements by id and connectivity.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("ELSET")
                    .optional("EALL")
                    .doc("Element set to activate while reading")
                .key("TYPE")
                    .required()
                    .doc("Element topology/type")
                    .allowed({
                        "C3D4", "C3D5", "C3D6", "C3D8", "C3D10", "C3D15", "C3D20", "C3D20R",
                        "B33", "T3", "P", "S3", "S4", "MITC4", "S6", "S8"
                    })
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& elset = keys.raw("ELSET");
            model._data->elem_sets.activate(elset);
        });

#define FEM_ADD_ELEMENT_VARIANT(KEY, TYPE, COUNT, DESC) \
        command.variant(fem::dsl::Variant::make() \
            .when(fem::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::dsl::Segment::make() \
                .range(fem::dsl::LineRange{}.min(1)) \
                .pattern(fem::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<fem::ID>().name("ID").desc("Element id") \
                    .fixed<fem::ID, COUNT>().name("N").desc(DESC) \
                ) \
                .bind([&model](fem::ID id, const std::array<fem::ID, COUNT>& nodes) { \
                    detail::set_element_from_array<TYPE, COUNT>(model, id, nodes); \
                }) \
            ) \
        );

        FEM_ADD_ELEMENT_VARIANT("C3D4", fem::model::C3D4, 4, "C3D4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D6", fem::model::C3D6, 6, "C3D6 connectivity (6 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D8", fem::model::C3D8, 8, "C3D8 connectivity (8 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D10", fem::model::C3D10, 10, "C3D10 connectivity (10 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D15", fem::model::C3D15, 15, "C3D15 connectivity (15 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D20", fem::model::C3D20, 20, "C3D20 connectivity (20 nodes)");
        FEM_ADD_ELEMENT_VARIANT("C3D20R", fem::model::C3D20R, 20, "C3D20R connectivity (20 nodes)");
        FEM_ADD_ELEMENT_VARIANT("B33", fem::model::B33, 2, "B33 connectivity (2 nodes)");
        FEM_ADD_ELEMENT_VARIANT("T3", fem::model::T3, 2, "T3 connectivity (2 nodes)");
        FEM_ADD_ELEMENT_VARIANT("P", fem::model::Point, 1, "Point element node");
        FEM_ADD_ELEMENT_VARIANT("S3", fem::model::S3, 3, "S3 connectivity (3 nodes)");
        FEM_ADD_ELEMENT_VARIANT("S4", fem::model::S4, 4, "S4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("MITC4", fem::model::MITC4, 4, "MITC4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_VARIANT("S6", fem::model::S6, 6, "S6 connectivity (6 nodes)");
        FEM_ADD_ELEMENT_VARIANT("S8", fem::model::S8, 8, "S8 connectivity (8 nodes)");

#undef FEM_ADD_ELEMENT_VARIANT

        // Special case: C3D5 mapped onto C3D8 with repeated apex node
        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"C3D5"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::ID>().name("ID").desc("Element id")
                    .fixed<fem::ID, 5>().name("N").desc("C3D5 connectivity (5 nodes)")
                )
                .bind([&model](fem::ID id, const std::array<fem::ID, 5>& nodes) {
                    model.set_element<fem::model::C3D8>(id,
                        nodes[0], nodes[1], nodes[2], nodes[3],
                        nodes[4], nodes[4], nodes[4], nodes[4]);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
