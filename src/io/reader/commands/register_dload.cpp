/**
 * @file register_dload.cpp
 * @brief Registers native and Abaqus-compatible distributed loads.
 *
 * FEMaster rows prescribe a traction vector directly on a compiled surface
 * region as TARGET, tx, ty, tz. Abaqus-compatible rows address elements and
 * dispatch by load label: BX/BY/BZ create volume force densities, GRAV creates
 * a density-scaled gravity load, and P<n> creates normal pressure. TRVEC<n>
 * and follower modes are currently rejected by the input reader.
 *
 * Named definitions select an explicit collector through Parser. Native
 * unnamed tractions are step-local; Abaqus-format DLOAD rows update source history.
 * Command entry resolves shared amplitudes. Nodal TRANSFORM assignments do not
 * affect distributed loads;
 * ORIENTATION applies only to traction vectors whose components require a local
 * basis.
 *
 * @see bc::DLoad
 * @see bc::PLoad
 * @see bc::VLoad
 * @see bc::InertialLoad
 * @see Parser::activate_load_collector
 * @see model::Model::resolve_element_region
 * @see model::Model::resolve_surface_region
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "../parser.h"
#include "../../../bc/structural/load_d.h"
#include "../../../bc/structural/load_inertial.h"
#include "../../../bc/structural/load_p.h"
#include "../../../bc/structural/load_v.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../../model/geometry/surface/surface_interface.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers surface tractions and element-based distributed loads.
 *
 * Native FEMaster rows resolve a surface target and store the supplied traction
 * components directly in bc::DLoad. Abaqus-compatible rows resolve an element
 * target and translate the supported distributed-load labels to the matching
 * FEMaster load implementation.
 *
 * BX/BY/BZ are global force densities and therefore map to bc::VLoad without an
 * orientation. GRAV maps the prescribed gravity vector to the translational
 * acceleration of bc::InertialLoad so density scaling remains part of the load
 * integration. An empty GRAV target includes all compiled elements and auxiliary
 * point masses; explicit targets select only their element region. P<n> and
 * P<n> extracts face n from every target element and creates a private compiled
 * surface region for bc::PLoad. TRVEC<n> is currently rejected.
 *
 * NAME identifies reusable collector storage and is required outside an analysis.
 * Unnamed analysis rows replace matching DLOAD targets; OP=NEW first clears only
 * that family. Complex-valued loading
 * is rejected. TRVEC is not currently supported, and pressure remains
 * prohibited in nonlinear analyses without load-stiffness support.
 *
 * @param registry Registry receiving the command grammar and callbacks.
 * @param parser Parser supplying the compiled model and active analysis.
 */
void register_dload(dsl::Registry& registry, Parser& parser) {
    registry.command("DLOAD", [&](dsl::Command& command) {
        // Distributed loads operate on compiled regions in model or analysis scope.
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc(
            "Create FEMaster surface tractions or Abaqus-compatible element-based "
            "BX/BY/BZ, GRAV and P<n> loads (TRVEC unsupported)."
        );

        // Resolved modifiers are reset on entry and shared by this occurrence's rows.
        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);
        auto amplitude   = std::make_shared<bc::Amplitude::Ptr        >(nullptr);

        // NAME is the FEMaster collector extension. The remaining options are
        // shared by native and supported Abaqus-compatible distributed loads.
        // An explicit NAME selects reusable definition storage. Unnamed analysis
        // rows update direct history and never inherit a prior insertion target.
        auto collector = std::make_shared<bc::LoadCollector::Ptr>(nullptr);

        command.keyword(
            dsl::KeywordSpec::make()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                .key("NAME").alternative("LOAD_COLLECTOR").alternative("LOADCOLLECTOR")
                    .optional().doc("Collector name; optional inside LOADCASE")
                .key("ORIENTATION").optional().doc("Coordinate system for traction components")
                .key("AMPLITUDE"  ).optional().doc("Amplitude scaling the complete load")
                .key("FOLLOWER"   ).optional("NO").allowed({"NO"})
                .flag("REAL")
                .flag("IMAGINARY")
        );

        command.on_enter([&parser, orientation, amplitude, collector](const dsl::Keys& keys) {
            auto& model = parser.model();

            orientation->reset();
            amplitude  ->reset();

            const std::string orientation_name = keys.raw("ORIENTATION");
            const std::string amplitude_name   = keys.raw("AMPLITUDE");
            const std::string collector_name   = keys.raw("NAME");

            const auto* loadcase = parser.active_loadcase();
            logging::error(loadcase == nullptr || loadcase->type_name() != "EIGENFREQ",
                "DLOAD: not supported in a FREQUENCY step");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "DLOAD: REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "DLOAD: IMAGINARY is not supported");
            logging::error(orientation_name.empty() || model._data->coordinate_systems.has(orientation_name),
                "DLOAD: coordinate system ", orientation_name, " does not exist");
            logging::error(amplitude_name.empty() || model._data->amplitudes.has(amplitude_name),
                "DLOAD: amplitude ", amplitude_name, " does not exist");

            logging::error(!collector_name.empty() || loadcase != nullptr,
                "DLOAD: definitions outside an analysis require NAME");

            if (!orientation_name.empty())
                *orientation = model._data->coordinate_systems.get(orientation_name);
            if (!amplitude_name.empty())
                *amplitude = model._data->amplitudes.get(amplitude_name);

            collector->reset();
            if (!collector_name.empty()) {
                parser.activate_load_collector(collector_name);
                *collector = model._data->load_cols.get();
            }
            if (collector_name.empty() && keys.raw("OP") == "NEW")
                parser.clear_conditions(bc::DLOAD);
        });

        // ---------------------------------------------------------------------
        // Resolve Abaqus element-face labels into private surface regions
        // ---------------------------------------------------------------------

        const auto parse_face_id = [](const std::string& type, const std::string& prefix) {
            logging::error(type.size() > prefix.size(),
                "DLOAD: load type ", type, " requires a one-based element-face number");

            ID face_id{};
            const char* begin = type.data() + prefix.size();
            const char* end   = type.data() + type.size();
            const auto [ptr, ec] = std::from_chars(begin, end, face_id);

            logging::error(ec == std::errc{} && ptr == end && face_id > 0,
                "DLOAD: invalid element-face load type ", type);
            return face_id;
        };

        const auto materialize_face_region = [&parser](const std::string& target, ID face_id) {
            auto& model    = parser.model();
            auto  elements = model.resolve_element_region(target);

            logging::error(elements != nullptr && elements->size() > 0,
                "DLOAD: target element region ", target, " is empty");

            // Build and validate every requested face before changing compiled
            // surface storage so a bad element cannot leave a partial region.
            std::vector<model::SurfacePtr> surfaces;
            surfaces.reserve(elements->size());

            for (const ID element_id : *elements) {
                logging::error(element_id >= 0 && element_id < static_cast<ID>(model._data->elements.size()),
                    "DLOAD: compiled element ", element_id, " is out of range");

                auto& element = model._data->elements[static_cast<std::size_t>(element_id)];
                logging::error(element != nullptr,
                    "DLOAD: compiled element ", element_id, " is not configured");

                auto surface = element->surface(face_id);
                logging::error(surface != nullptr,
                    "DLOAD: face ", face_id, " is not available for element ", element_id, " (", element->type_name(), ")");

                surfaces.push_back(std::move(surface));
            }

            // Reuse matching compiled faces so repeated P<n> input
            // resolves to stable surface identifiers. Without this geometric
            // reuse MOD would append loads on fresh IDs instead of replacing
            // the same physical face in condition history.
            auto region = std::make_shared<model::SurfaceRegion>("INTERNAL");
            for (auto& surface : surfaces) {
                const auto existing = std::find_if(
                    model._data->surfaces.begin(),
                    model._data->surfaces.end(),
                    [&](const model::SurfacePtr& candidate) {
                        return candidate != nullptr && candidate->n_nodes == surface->n_nodes
                            && std::equal(surface->begin(), surface->end(), candidate->begin());
                    }
                );

                // Connectivity ordering includes the face orientation and thus
                // preserves the pressure-normal and traction integration basis.
                ID surface_id = static_cast<ID>(std::distance(model._data->surfaces.begin(), existing));
                if (existing == model._data->surfaces.end()) {
                    model._data->surfaces.push_back(std::move(surface));
                }
                region->add(surface_id);
            }

            return region;
        };

        // ---------------------------------------------------------------------
        // FEMaster format: SURFACE, tx, ty, tz
        // ---------------------------------------------------------------------

        command.data(
            dsl::Pattern::make()
                .one<std::string   >("TARGET", "Compiled surface set or scalar reference")
                .fixed<Precision, 3>("LOAD", "Traction components tx, ty, tz")
                    .defaults(Precision{0}),
            [&parser, orientation, amplitude, collector](
                const std::string&              target,
                const std::array<Precision, 3>& values
            ) {
                auto& model = parser.model();

                // Native DLOAD addresses an existing surface region directly.
                auto load = std::make_shared<bc::DLoad>();
                load->region_      = model.resolve_surface_region(target);
                load->values_      = Vec3{values[0], values[1], values[2]};
                load->orientation_ = *orientation;
                load->amplitude_   = *amplitude;
                // Named definitions and direct history are mutually exclusive targets
                if (*collector) {
                    (*collector)->add(std::move(load));
                } else {
                    parser.select_collector_condition(bc::DLOAD, std::move(load));
                }
            }
        );

        // ---------------------------------------------------------------------
        // Abaqus format: ELEMENT, TYPE, MAGNITUDE [, dx, dy, dz]
        // ---------------------------------------------------------------------

        command.data(
            dsl::Pattern::make()
                .one<std::string   >("TARGET", "Compiled element set or scalar element reference")
                    .on_empty  (std::string{})
                .one<std::string   >("TYPE", "Abaqus distributed-load label")
                .one<Precision     >("MAGNITUDE", "Load magnitude")
                .fixed<Precision, 3>("DIRECTION", "Optional direction components")
                    .defaults(std::numeric_limits<Precision>::quiet_NaN()),
            [&parser, orientation, amplitude, parse_face_id, materialize_face_region, collector](
                const std::string&              target,
                const std::string&              type,
                Precision                       magnitude,
                const std::array<Precision, 3>& direction
            ) {
                auto& model = parser.model();

                // Preserve whether TARGET was omitted. All-element GRAV also
                // includes auxiliary POINTMASS objects outside the ELSET namespace;
                // an explicit target, including EALL, selects only its region.
                const bool        all_elements   = target.empty();
                const std::string element_target = all_elements ? "EALL" : target;
                const std::string identifier = "ELEMENT:" + element_target + ":" + type;

                const bool direction_omitted  =  std::isnan(direction[0]) &&  std::isnan(direction[1]) &&  std::isnan(direction[2]);
                const bool direction_complete = !std::isnan(direction[0]) && !std::isnan(direction[1]) && !std::isnan(direction[2]);

                logging::error(type.compare(0, 5, "TRVEC") != 0,
                    "DLOAD: TRVEC is not supported");
                logging::error(std::isfinite(magnitude),
                    "DLOAD: magnitude must be finite");

                // BX/BY/BZ are force densities in the fixed global basis.
                // They map directly to the corresponding component of VLoad.
                if (type == "BX" || type == "BY" || type == "BZ") {
                    logging::error(direction_omitted,
                        "DLOAD: ", type, " accepts no direction components");

                    // Preserve the DLOAD component as history identity.
                    // NaN marks axes not prescribed by this BX/BY/BZ entry;
                    // VLoad converts them to zero only during assembly.
                    Vec3 body_force = Vec3::Constant(NAN);
                    if (type == "BX") body_force[0] = magnitude;
                    if (type == "BY") body_force[1] = magnitude;
                    if (type == "BZ") body_force[2] = magnitude;

                    auto load = std::make_shared<bc::VLoad>();
                    load->region_    = model.resolve_element_region(element_target);
                    load->values_    = body_force;
                    load->amplitude_ = *amplitude;
                    // Named definitions and direct history are mutually exclusive targets
                    if (*collector) {
                        (*collector)->add(std::move(load));
                    } else {
                        parser.modify_conditions(bc::DLOAD, identifier, {std::move(load)});
                    }
                    return;
                }

                // GRAV is an acceleration and therefore uses density-scaled
                // inertia integration. InertialLoad applies -rho*a, so the
                // prescribed gravity vector enters as the negative body
                // acceleration to reproduce +rho*g as the external load.
                if (type == "GRAV") {
                    const Vec3 gravity_direction{direction[0], direction[1], direction[2]};
                    const Precision direction_norm = gravity_direction.norm();

                    logging::error(direction_complete,
                        "DLOAD: GRAV requires three direction components");
                    logging::error(direction_norm > Precision(0) && std::isfinite(direction_norm),
                        "DLOAD: GRAV direction must be finite and nonzero");

                    auto load = std::make_shared<bc::InertialLoad>();
                    load->region_                = model.resolve_element_region(element_target);
                    load->center_                = Vec3::Zero();
                    load->center_acc_            = -magnitude * gravity_direction;
                    load->omega_                 = Vec3::Zero();
                    load->alpha_                 = Vec3::Zero();
                    load->amplitude_             = *amplitude;
                    load->consider_point_masses_ = all_elements;
                    // Named definitions and direct history are mutually exclusive targets
                    if (*collector) {
                        (*collector)->add(std::move(load));
                    } else {
                        parser.modify_conditions(bc::DLOAD, identifier, {std::move(load)});
                    }
                    return;
                }

                const auto* loadcase  = parser.active_loadcase();
                const bool nonlinear = loadcase != nullptr && loadcase->type_name() == "NONLINEARSTATIC";

                // P<n> addresses face n of every target element. Pressure is
                // follower loading in Abaqus; the current PLoad does not
                // provide the nonlinear load-stiffness contribution.
                if (type.size() > 1 && type[0] == 'P') {
                    logging::error(direction_omitted,
                        "DLOAD: ", type, " accepts no direction components");
                    logging::error(!nonlinear,
                        "DLOAD: follower pressure is not supported in nonlinear steps");

                    const ID face_id = parse_face_id(type, "P");

                    auto load = std::make_shared<bc::PLoad>();
                    load->region_    = materialize_face_region(element_target, face_id);
                    load->pressure_  = magnitude;
                    load->amplitude_ = *amplitude;
                    // Named definitions and direct history are mutually exclusive targets
                    if (*collector) {
                        (*collector)->add(std::move(load));
                    } else {
                        parser.modify_conditions(bc::DLOAD, identifier, {std::move(load)});
                    }
                    return;
                }

                logging::error(false,
                    "DLOAD: supported Abaqus load types are BX, BY, BZ, GRAV and P<n>");
            }
        );
    });
}

} // namespace fem::io::reader::commands
