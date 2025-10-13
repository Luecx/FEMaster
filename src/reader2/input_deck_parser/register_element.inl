// commands/register_element.inl
// Intentionally no include guards — include this in exactly ONE translation unit.

// --- element type headers (adjust paths to your project layout) ---
#include "../../model/beam/b33.h"
#include "../../model/truss/truss.h"
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

#include <tuple>   // std::apply

namespace fem::reader2::commands {

inline void register_element(Registry& reg, fem::model::Model& model, reader::Writer& writer) {
    (void)writer; // reserved for future use

    reg.register_command("ROOT", "ELEMENT",
        [&model](CommandSpec& c) {
            c.doc("*ELEMENT: Define finite elements by ID and connectivity. "
                  "Connectivity may span multiple data lines.");

            // keyword args
            c.keys(KeyRules::make()
                .optional("ELSET", "EALL").describe("ELSET", "Active element set")
                .require("TYPE").allowed("TYPE", {
                    // solids
                    "C3D4","C3D5","C3D6","C3D8","C3D10","C3D15","C3D20","C3D20R",
                    // beams/truss/point
                    "B33","T3","P",
                    // shells
                    "S3","S4","MITC4","S6","S8"
                }).describe("TYPE","Element topology/type")
            );

            // per-keyword hook: activate ELSET
            c.on_keyword([&model](Context&, const Keyword& kw) {
                model._data->elem_sets.activate(kw.get<std::string>("ELSET"));
            });

            // ------------------------------------------------------------------
            // X-macro list of (TYPE string, Element class, Node count)
            // ------------------------------------------------------------------
#define FEM_ELEM_LIST \
    X("C3D4",   fem::model::C3D4,   4)  \
    X("C3D6",   fem::model::C3D6,   6)  \
    X("C3D8",   fem::model::C3D8,   8)  \
    X("C3D10",  fem::model::C3D10, 10)  \
    X("C3D15",  fem::model::C3D15, 15)  \
    X("C3D20",  fem::model::C3D20, 20)  \
    X("C3D20R", fem::model::C3D20R, 20) \
    X("B33",    fem::model::B33,    2)  \
    X("T3",     fem::model::T3,     2)  \
    X("P",      fem::model::Point,  1)  \
    X("S3",     fem::model::S3,     3)  \
    X("S4",     fem::model::S4,     4)  \
    X("MITC4",  fem::model::MITC4,  4)  \
    X("S6",     fem::model::S6,     6)  \
    X("S8",     fem::model::S8,     8)

            // ------------------------------------------------------------------
            // Generic generator for all "normal" types in FEM_ELEM_LIST
            // (no outside helpers; use std::apply directly)
            // ------------------------------------------------------------------
#define X(TYPESTR, ELEMCLS, NN) \
    c.plan( \
        Condition::eq("TYPE", TYPESTR), \
        Segment::make() \
            .range(Range::make().min(1)) \
            .pattern( \
                Pattern::make() \
                    .allow_multiline(true) \
                    .fixed<int,1>("ID").desc("element id") \
                    .fixed<int,NN>("N").desc("connectivity (" TYPESTR " : " #NN " nodes)") \
                    .bind_kv<int, std::array<int, NN>>( \
                        [&model](Context&, const Keyword::Map&, int id, const std::array<int, NN>& n, const LineMeta&) { \
                            std::apply([&](auto... nn){ \
                                model.set_element<ELEMCLS>(id, nn...); \
                            }, n); \
                        } \
                    ) \
            ) \
            .label(TYPESTR) \
    );

            FEM_ELEM_LIST
#undef X
#undef FEM_ELEM_LIST

            // ------------------------------------------------------------------
            // Special-case: C3D5 is mapped to C3D8 by repeating the 5th node
            // ------------------------------------------------------------------
            c.plan(
                Condition::eq("TYPE","C3D5"),
                Segment::make()
                    .range(Range::make().min(1))
                    .pattern(
                        Pattern::make()
                            .allow_multiline(true)
                            .fixed<int,1>("ID").desc("element id")
                            .fixed<int,5>("N").desc("connectivity (5 nodes, mapped to C3D8)")
                            .bind_kv<int, std::array<int,5>>(
                                [&model](Context&, const Keyword::Map&, int id, const std::array<int,5>& n, const LineMeta&) {
                                    model.set_element<fem::model::C3D8>(
                                        id,
                                        n[0], n[1], n[2], n[3],
                                        n[4], n[4], n[4], n[4]
                                    );
                                })
                    )
                    .label("C3D5 (mapped to C3D8)")
            );
        }
    );
}

} // namespace fem::reader2::commands
