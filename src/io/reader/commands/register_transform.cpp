/**
 * @file register_transform.cpp
 * @brief Registers nodal coordinate transformations.
 *
 * `*TRANSFORM` executes after `Model::compile()` so the referenced node set
 * already contains dense compiled node identifiers. It may appear at root or
 * assembly scope; retaining the assembly parent prevents a transform from
 * prematurely closing a still-active assembly during the post-compile replay.
 *
 * Rectangular (`TYPE=R`) and cylindrical (`TYPE=C`) transformations are
 * supported. Coordinate-system construction is parser-owned; `Model` only
 * registers the finished polymorphic object.
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <limits>
#include <memory>
#include <string>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../cos/cylindrical_system.h"
#include "../../../cos/rectangular_system.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

/**
 * @brief Registers the *TRANSFORM command.
 *
 * The command assigns a rectangular or cylindrical local coordinate system to
 * all nodes of a compiled node set. The nodal degrees of freedom themselves
 * remain in the global coordinate system; the resolved node-to-coordinate-system
 * mapping is stored in the parser and later applied when loads and supports are
 * materialized.
 *
 * Rectangular transformations use the two supplied points to define the local
 * basis. Cylindrical transformations interpret the points as the cylinder axis
 * and construct the corresponding FEMaster cylindrical coordinate system.
 *
 * @param registry Command registry receiving the TRANSFORM grammar.
 * @param parser   Parser storing the resolved nodal transformation assignments.
 */
void register_transform(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("TRANSFORM", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Assign a rectangular or cylindrical transform to a compiled node set.");

        auto nset = std::make_shared<std::string>();
        auto type = std::make_shared<std::string>();

        // ---------------------------------------------------------------------
        // Parse the transform target and transformation type
        // ---------------------------------------------------------------------

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET").required().doc("Compiled node set receiving the transformation")
                .key("TYPE").optional("R").allowed({"R", "C"}).doc("R = rectangular, C = cylindrical")
        );

        command.on_enter([nset, type](const fem::io::dsl::Keys& keys) {
            *nset = keys.raw("NSET");
            *type = keys.raw("TYPE");
        });

        // ---------------------------------------------------------------------
        // Parse the two points defining the local coordinate system
        // ---------------------------------------------------------------------

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 6>().name("DATA").desc("Coordinates of points a and b")
                )
                .bind([&parser, nset, type](const std::array<fem::Precision, 6>& data) {
                    auto& model      = parser.model();
                    auto& transforms = parser.node_transforms;

                    const fem::Vec3      a{data[0], data[1], data[2]};
                    const fem::Vec3      b{data[3], data[4], data[5]};
                    const fem::Precision eps = std::numeric_limits<fem::Precision>::epsilon();
                    const std::string    orientation = "__TRANSFORM_" + *nset;

                    // Resolve the compiled node set before constructing the
                    // coordinate system assigned to its nodes.
                    logging::error(model._data->node_sets.has(*nset),
                        "TRANSFORM: node set ", *nset, " is not defined");
                    logging::error(!model._data->coordinate_systems.has(orientation),
                        "TRANSFORM: node set ", *nset, " is defined more than once");

                    // Construct the FEMaster coordinate-system representation
                    // corresponding to the requested TRANSFORM type.
                    if (*type == "R") {
                        const fem::Precision norm_a = a.norm();
                        const fem::Precision norm_b = b.norm();
                        const fem::Precision cross  = a.cross(b).norm();

                        logging::error(norm_a > eps && norm_b > eps && cross > eps * norm_a * norm_b,
                            "TRANSFORM TYPE=R requires nonzero, non-collinear points a and b");

                        model.add_coordinate_system(
                            std::make_shared<cos::RectangularSystem>(orientation, a, b)
                        );
                    } else {
                        const fem::Vec3 axis = b - a;

                        logging::error(axis.norm() > eps,
                            "TRANSFORM TYPE=C requires distinct axis points a and b");

                        const fem::Vec3 axial      = axis.normalized();
                        const fem::Vec3 radial     = axial.unitOrthogonal();
                        const fem::Vec3 tangential = axial.cross(radial).normalized();

                        model.add_coordinate_system(
                            std::make_shared<cos::CylindricalSystem>(
                                orientation, a, a + radial, a + tangential)
                        );
                    }

                    // Store the resolved nodal assignment. Loads and supports
                    // later use this mapping instead of transforming the nodal
                    // degrees of freedom themselves.
                    auto nodes = model._data->node_sets.get(*nset);

                    for (const fem::ID node_id : *nodes) {
                        auto [it, inserted] = transforms.emplace(node_id, orientation);

                        logging::error(inserted || it->second == orientation,
                            "TRANSFORM: node ", node_id, " belongs to multiple incompatible definitions");
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands