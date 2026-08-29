/**
 * @file register_point_mass.cpp
 * @brief Registers the native FEMaster `POINTMASS` command on compiled node sets.
 *
 * The native FEMaster point-mass syntax is node-set based: `POINTMASS, NSET=...`
 * assigns concentrated translational mass, rotary inertia and optional diagonal
 * translational/rotational ground stiffness to every node of an existing NSET.
 * The command is evaluated after `Model::compile()`, so the referenced NSET and
 * all node ids already belong to the dense compiled assembly namespace.
 *
 * Point-mass physics is represented through the same element/section mechanism
 * used by Abaqus MASS elements. For every target node this command creates one
 * auxiliary `model::PointElement`; all generated elements share one
 * `PointMassSection` containing the values supplied by the command. The former
 * non-element PointMass feature is therefore not involved in this path.
 *
 * These auxiliary point elements are created after the dense model topology has
 * been frozen. They are intentionally stored in `ModelData::point_elements`
 * rather than inserted into `ModelData::elements` or any compiled ELSET. This
 * preserves the compiled element numbering and element-domain field offsets while
 * still allowing the point elements to participate in active-DOF, stiffness,
 * mass and inertia assembly through the normal `PointElement` implementation.
 *
 * The public command syntax remains NSET-based; the generated PointElements and
 * their internal ElementRegion are implementation details used to express the
 * nodal property through the common section architecture.
 *
 * @see model::PointElement
 * @see PointMassSection
 * @see model::ModelData::point_elements
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <memory>
#include <string>

#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../model/element/point.h"
#include "../../../model/model.h"
#include "../../../section/section_point_mass.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers `POINTMASS, NSET=...` in the native FEMaster reader.
 *
 * The command is valid only at ROOT scope and operates on the compiled model.
 * Its data record contains one scalar translational mass, three rotary inertia
 * components, three translational ground-spring constants and three rotational
 * ground-spring constants. Missing values retain the existing zero defaults.
 *
 * The bind resolves the compiled target NSET, constructs one shared
 * `PointMassSection`, then creates one auxiliary PointElement for every node in
 * the set. Each generated element references exactly one compiled node and the
 * shared section. Negative internal element ids are used only to distinguish
 * these post-compile objects from each other; they are not inserted into the
 * compiled element namespace or exposed through ELSETs.
 *
 * @param registry Parser registry receiving the native POINTMASS command.
 * @param model Compiled model receiving the auxiliary point elements.
 */
void register_point_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("POINTMASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign concentrated mass, inertia and ground stiffness to a compiled node set.");

        // Keep the target NSET name as command-local state. The keyword line is
        // parsed before the data variant, so the bind later resolves this name
        // against the already compiled node-set namespace.
        auto nset = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET").required().doc("Target compiled node set")
        );

        command.on_enter([nset](const fem::io::dsl::Keys& keys) {
            *nset = keys.raw("NSET");
        });

        // ROOT variant -- create auxiliary PointElements on the compiled model.
        //
        // Native POINTMASS has only one semantic scope. `parent_is("ROOT")`
        // expresses that restriction directly in the DSL, so the bind does not
        // need to inspect any parser scope state. The parser reaches this command
        // after Model::compile(), which means NSET members are already dense
        // assembly node ids. The callback validates that compiled state, converts
        // the four command data groups into the corresponding section values,
        // creates one shared PointMassSection, and finally creates one one-node
        // PointElement per target node. The generated elements are retained in
        // ModelData::point_elements instead of the frozen dense element array.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::parent_is("ROOT"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::Precision>().name("MASS").desc("Point mass")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("INERTIA").desc("Rotary inertia Ix,Iy,Iz")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("SPRING").desc("Translational spring constants")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("ROTSPRING").desc("Rotational spring constants")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, nset](
                          fem::Precision mass,
                          const std::array<fem::Precision, 3>& inertia_data,
                          const std::array<fem::Precision, 3>& spring_data,
                          const std::array<fem::Precision, 3>& rotary_data) {
                    // Resolve the public NSET target in compiled assembly space.
                    // POINTMASS is intentionally a post-compile command and must
                    // never create semantic Part-local topology.
                    logging::error(model._data != nullptr && model._data->compiled,
                        "POINTMASS: command requires a compiled model");
                    logging::error(model._data->node_sets.has(*nset),
                        "POINTMASS: node set ", *nset, " does not exist");

                    const auto nodes = model._data->node_sets.get(*nset);
                    logging::error(nodes != nullptr,
                        "POINTMASS: node set ", *nset, " is not initialized");

                    // Convert the parsed fixed-size arrays into the Vec3 values
                    // stored by PointMassSection. Each vector acts diagonally in
                    // the compiled global nodal directions.
                    const fem::Vec3 inertia{inertia_data[0], inertia_data[1], inertia_data[2]};
                    const fem::Vec3 spring {spring_data[0],  spring_data[1],  spring_data[2]};
                    const fem::Vec3 rotary {rotary_data[0],  rotary_data[1],  rotary_data[2]};

                    // All PointElements created by this command share one section
                    // because the POINTMASS record assigns identical properties to
                    // every node in the target NSET. The internal region records
                    // those auxiliary element ids only for section ownership; it
                    // is not registered as a public compiled ELSET.
                    auto section_region = std::make_shared<model::ElementRegion>(
                        "__POINTMASS_" + std::to_string(model._data->point_elements.size()));
                    auto section = std::make_shared<PointMassSection>(
                        section_region, mass, inertia, spring, rotary);

                    // Materialize one independent one-node element per NSET node.
                    // Negative ids keep these post-compile elements distinct
                    // without colliding with the non-negative dense element ids
                    // established by compile(). They remain outside all element-
                    // domain enumeration and are assembled through point_elements.
                    for (const ID node : *nodes) {
                        const ID id = -1 - static_cast<ID>(model._data->point_elements.size());
                        section_region->add(id);

                        auto point = std::make_shared<model::PointElement>(
                            id, std::array<ID, model::PointElement::N>{node});
                        point->_model_data = model._data.get();
                        point->set_section(section);
                        model._data->point_elements.push_back(std::move(point));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
