/**
 * @file model.cpp
 * @brief Implements high-level finite-element model construction operations.
 *
 * This file implements the `Model` facade operations that connect named regions,
 * coordinate systems, loads, supports, constraints, sections and features with
 * the persistent repositories owned by `ModelData`. Element-local mechanics and
 * global matrix/field assembly remain in the corresponding element classes and
 * `model_build.cpp`.
 *
 * Surface-to-surface contact construction is intentionally narrow: both MASTER
 * and SLAVE names resolve to non-empty surface regions and are forwarded with the
 * user penalty, clearance and normal orientation to `constraint::Contact`.
 * Geometric mortar processing itself remains entirely inside the contact class.
 *
 * @see Model
 * @see ModelData
 * @see constraint::Contact
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#include "model.h"

#include "../bc/load_collector.h"
#include "../bc/support_collector.h"
#include "../feature/point_mass.h"
#include "../section/section.h"
#include "../section/section_beam.h"
#include "../section/section_shell_abd.h"
#include "../section/section_shell_integrated.h"
#include "../section/section_solid.h"
#include "../section/section_truss.h"

namespace fem {
namespace model {

/**
 * Initializes step-local state and caches of all structural elements.
 *
 * If any element initialization throws, already initialized elements are cleaned
 * up through `step_end()` before the exception is propagated. Non-structural and
 * empty element slots are ignored.
 */
void Model::step_begin() {
    try {
        for (auto& elem : _data->elements) {
            if (elem != nullptr) {
                if (auto* structural = elem->as<StructuralElement>()) {
                    structural->step_begin();
                }
            }
        }
    } catch (...) {
        step_end();
        throw;
    }
}

/**
 * Releases step-local state and caches of all structural elements.
 *
 * The operation mirrors `step_begin()` and ignores empty/non-structural entries.
 */
void Model::step_end() {
    for (auto& elem : _data->elements) {
        if (elem != nullptr) {
            if (auto* structural = elem->as<StructuralElement>()) {
                structural->step_end();
            }
        }
    }
}

/**
 * Adds a two-node connector defined by singleton node sets and one coordinate
 * system. All referenced entities are validated before the connector is stored.
 */
void Model::add_connector(const std::string& set1,
                          const std::string& set2,
                          const std::string& coordinate_system,
                          constraint::ConnectorType type) {
    logging::error(_data->node_sets.has(set1),
        "Node set ", set1, " does not exist");
    logging::error(_data->node_sets.has(set2),
        "Node set ", set2, " does not exist");
    logging::error(_data->coordinate_systems.has(coordinate_system),
        "Coordinate system ", coordinate_system, " does not exist");
    logging::error(_data->node_sets.get(set1)->size() == 1,
        "Set 1 must contain exactly one node");
    logging::error(_data->node_sets.get(set2)->size() == 1,
        "Set 2 must contain exactly one node");

    const ID id1 = _data->node_sets.get(set1)->first();
    const ID id2 = _data->node_sets.get(set2)->first();

    _data->connectors.emplace_back(
        id1,
        id2,
        _data->coordinate_systems.get(coordinate_system),
        type
    );
}

/**
 * Adds a kinematic or distributing coupling from one master node to a node or
 * surface slave region.
 *
 * The master region must be a singleton node set. `is_surface` selects whether
 * `slave_set` is resolved from the surface or node-set repository. The master
 * region pointer is retained on the created coupling for later reporting.
 */
void Model::add_coupling(const std::string& master_set,
                         const std::string& slave_set,
                         Dofs coupled_dofs,
                         constraint::CouplingType type,
                         bool is_surface) {
    logging::error(_data->node_sets.get(master_set)->size() == 1,
        "Master set must contain exactly one node");

    if (is_surface) {
        logging::error(_data->surface_sets.has(slave_set),
            "Slave set ", slave_set, " is not a defined surface set");
    } else {
        logging::error(_data->node_sets.has(slave_set),
            "Slave set ", slave_set, " is not a defined node set");
    }

    const ID master_node = _data->node_sets.get(master_set)->first();

    if (is_surface) {
        _data->couplings.emplace_back(
            master_node,
            _data->surface_sets.get(slave_set),
            coupled_dofs,
            type
        );
    } else {
        _data->couplings.emplace_back(
            master_node,
            _data->node_sets.get(slave_set),
            coupled_dofs,
            type
        );
    }

    if (!_data->couplings.empty()) {
        auto& coupling = _data->couplings.back();
        coupling.master_region = _data->node_sets.get(master_set);
    }
}

/**
 * Adds a tie constraint between a master surface/line region and a node/surface
 * slave region.
 *
 * Slave node sets are preferred when the same name exists in both repositories;
 * otherwise a non-empty surface set is used. The master resolves first to a
 * non-empty surface region and then to a non-empty line region. The geometric
 * tie distance and optional adjustment flag are forwarded unchanged.
 */
void Model::add_tie(const std::string& master_set,
                    const std::string& slave_set,
                    Precision distance,
                    bool adjust) {
    // ---------------------------------------------------------------------
    // Resolve slave: prefer node set, otherwise surface set
    // ---------------------------------------------------------------------
    NodeRegion::Ptr    slave_node_ptr    = nullptr;
    SurfaceRegion::Ptr slave_surface_ptr = nullptr;

    if (_data->node_sets.has(slave_set) && _data->node_sets.get(slave_set) &&
        _data->node_sets.get(slave_set)->size() > 0) {
        slave_node_ptr = _data->node_sets.get(slave_set);
    } else if (_data->surface_sets.has(slave_set) && _data->surface_sets.get(slave_set) &&
               _data->surface_sets.get(slave_set)->size() > 0) {
        slave_surface_ptr = _data->surface_sets.get(slave_set);
    } else {
        logging::error(false,
            "Slave set ", slave_set,
            " is neither a defined non-empty node set nor a defined non-empty surface set");
    }

    // ---------------------------------------------------------------------
    // Resolve master: surface region first, then line region
    // ---------------------------------------------------------------------
    const bool has_surfaces =
        _data->surface_sets.has(master_set) && _data->surface_sets.get(master_set) &&
        _data->surface_sets.get(master_set)->size() > 0;

    const bool has_lines =
        _data->line_sets.has(master_set) && _data->line_sets.get(master_set) &&
        _data->line_sets.get(master_set)->size() > 0;

    if (has_surfaces) {
        SurfaceRegion::Ptr master_ptr = _data->surface_sets.get(master_set);

        if (slave_node_ptr) {
            _data->ties.emplace_back(master_ptr, slave_node_ptr, distance, adjust);
        } else {
            _data->ties.emplace_back(master_ptr, slave_surface_ptr, distance, adjust);
        }
        return;
    }

    if (has_lines) {
        LineRegion::Ptr master_line_ptr = _data->line_sets.get(master_set);

        if (slave_node_ptr) {
            _data->ties.emplace_back(master_line_ptr, slave_node_ptr, distance, adjust);
        } else {
            _data->ties.emplace_back(master_line_ptr, slave_surface_ptr, distance, adjust);
        }
        return;
    }

    logging::error(false,
        "Master set ", master_set, " contains neither surfaces nor lines");
}

/**
 * Adds one frictionless dual-mortar surface-to-surface contact definition.
 *
 * Both region names must resolve to non-empty surface sets. The input `penalty`
 * is the starting augmented-Lagrange penalty; later stagnation-based adaptation
 * is numerical contact state and remains internal to `constraint::Contact`.
 * `clearance` is subtracted from the signed normal gap and `flip_normal` reverses
 * the master normal orientation during contact evaluation.
 */
void Model::add_contact(const std::string& master_set,
                        const std::string& slave_set,
                        Precision penalty,
                        Precision clearance,
                        bool flip_normal) {
    logging::error(_data->surface_sets.has(master_set)
                && _data->surface_sets.get(master_set)
                && _data->surface_sets.get(master_set)->size() > 0,
        "CONTACT: master set ", master_set,
        " is not a defined non-empty surface set");

    logging::error(_data->surface_sets.has(slave_set)
                && _data->surface_sets.get(slave_set)
                && _data->surface_sets.get(slave_set)->size() > 0,
        "CONTACT: slave set ", slave_set,
        " is not a defined non-empty surface set");

    _data->contacts.emplace_back(
        _data->surface_sets.get(master_set),
        _data->surface_sets.get(slave_set),
        penalty,
        clearance,
        flip_normal
    );
}

/**
 * Adds a rigid-body-motion suppression constraint for a non-empty element set.
 */
void Model::add_rbm(const std::string& set) {
    logging::error(_data->elem_sets.has(set),
        "RBM: element set ", set, " not found");
    logging::error(_data->elem_sets.get(set) && _data->elem_sets.get(set)->size() > 0,
        "RBM: element set ", set, " is empty");

    _data->rbms.emplace_back(_data->elem_sets.get(set));
}

/**
 * Adds a concentrated six-component load to a named node region.
 *
 * Optional coordinate-system and amplitude names are resolved before forwarding
 * the load to the currently active load collector.
 */
void Model::add_cload(const std::string& nset,
                      Vec6 load,
                      const std::string& orientation,
                      const std::string& amplitude) {
    logging::error(_data->node_sets.has(nset),
        "Node set ", nset, " does not exist");

    auto region_ptr = _data->node_sets.get(nset);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_cload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

/**
 * Adds a concentrated six-component load to one global node id.
 *
 * A temporary one-node region is created so the same collector representation is
 * used as for named node-set loads. Optional orientation and amplitude handling
 * is identical to the named-region overload.
 */
void Model::add_cload(ID id,
                      Vec6 load,
                      const std::string& orientation,
                      const std::string& amplitude) {
    auto region_ptr = std::make_shared<NodeRegion>("INTERNAL");
    region_ptr->add(id);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_cload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

/**
 * Adds a distributed vector surface load to a named surface region.
 */
void Model::add_dload(const std::string& sfset,
                      Vec3 load,
                      const std::string& orientation,
                      const std::string& amplitude) {
    logging::error(_data->surface_sets.has(sfset),
        "Surface set ", sfset, " does not exist");

    auto region_ptr = _data->surface_sets.get(sfset);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_dload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

/**
 * Adds a distributed vector surface load to one global surface id by constructing
 * a temporary one-surface region.
 */
void Model::add_dload(ID id,
                      Vec3 load,
                      const std::string& orientation,
                      const std::string& amplitude) {
    auto region_ptr = std::make_shared<SurfaceRegion>("INTERNAL");
    region_ptr->add(id);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_dload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

/**
 * Adds a scalar pressure load to a named surface region with an optional
 * amplitude.
 */
void Model::add_pload(const std::string& sfset,
                      Precision load,
                      const std::string& amplitude) {
    logging::error(_data->surface_sets.has(sfset),
        "Surface set ", sfset, " does not exist");

    auto region_ptr = _data->surface_sets.get(sfset);

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_pload(region_ptr, load, amplitude_ptr);
}

/**
 * Adds a scalar pressure load to one global surface id using a temporary surface
 * region.
 */
void Model::add_pload(ID id, Precision load, const std::string& amplitude) {
    auto region_ptr = std::make_shared<SurfaceRegion>("INTERNAL");
    region_ptr->add(id);

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_pload(region_ptr, load, amplitude_ptr);
}

/**
 * Adds a distributed vector volume load to a named element region.
 */
void Model::add_vload(const std::string& elset,
                      Vec3 load,
                      const std::string& orientation,
                      const std::string& amplitude) {
    logging::error(_data->elem_sets.has(elset),
        "Element set ", elset, " does not exist");

    auto region_ptr = _data->elem_sets.get(elset);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_vload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

/**
 * Adds a distributed vector volume load to one global element id using a
 * temporary element region.
 */
void Model::add_vload(ID id,
                      Vec3 load,
                      const std::string& orientation,
                      const std::string& amplitude) {
    auto region_ptr = std::make_shared<ElementRegion>("INTERNAL");
    region_ptr->add(id);

    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
        orientation_ptr = _data->coordinate_systems.get(orientation);
    }

    bc::Amplitude::Ptr amplitude_ptr = nullptr;
    if (!amplitude.empty()) {
        logging::error(_data->amplitudes.has(amplitude),
            "Amplitude ", amplitude, " does not exist");
        amplitude_ptr = _data->amplitudes.get(amplitude);
    }

    _data->load_cols.get()->add_vload(region_ptr, load, orientation_ptr, amplitude_ptr);
}

/**
 * Adds a rigid-body inertia load field to a named element region.
 *
 * The center, translational acceleration, angular velocity and angular
 * acceleration are forwarded to the active load collector. Point-mass features
 * participate only when explicitly requested.
 */
void Model::add_inertialload(const std::string& elset,
                             Vec3 center,
                             Vec3 center_acceleration,
                             Vec3 angular_velocity,
                             Vec3 angular_acceleration,
                             bool consider_point_masses) {
    logging::error(_data->elem_sets.has(elset),
        "Element set ", elset, " does not exist");

    auto region_ptr = _data->elem_sets.get(elset);
    _data->load_cols.get()->add_inertialload(
        region_ptr,
        center,
        center_acceleration,
        angular_velocity,
        angular_acceleration,
        consider_point_masses
    );
}

/**
 * Defines or redefines a named load amplitude and clears any previous samples.
 */
void Model::define_amplitude(const std::string& name, bc::Interpolation interpolation) {
    auto amplitude = _data->amplitudes.activate(name);
    amplitude->set_interpolation(interpolation);
    amplitude->clear_samples();
}

/**
 * Appends one `(time, value)` sample to an existing named amplitude.
 */
void Model::add_amplitude_sample(const std::string& name, Precision time, Precision value) {
    logging::error(_data->amplitudes.has(name),
        "Amplitude ", name, " has not been defined");

    _data->amplitudes.get(name)->add_sample(time, value);
}

/**
 * Adds a temperature load from a scalar nodal field and reference temperature.
 *
 * The temperature field must exist, use NODE domain and contain exactly one
 * component. Constitutive application remains the responsibility of the active
 * load collector and structural elements.
 */
void Model::add_tload(std::string& temp_field, Precision ref_temp) {
    auto field = _data->get_field(temp_field);

    logging::error(field != nullptr,
        "Temperature field ", temp_field, " does not exist");
    logging::error(field->domain == FieldDomain::NODE,
        "Temperature field ", temp_field, " must be a node field");
    logging::error(field->components == 1,
        "Temperature field ", temp_field, " must have 1 component");

    _data->load_cols.get()->add_tload(field, ref_temp);

    // // TODO
    // logging::error(_fields_temperature.has(temp_field), "Temperature field ", temp_field, " does not exist");
    //
    // auto temp_ptr = _fields_temperature.get(temp_field);
    // for (ElementPtr& elem : _data->elements) {
    //     if (elem == nullptr) continue;
    //     if (auto sel = elem->as<StructuralElement>())
    //         sel->apply_tload((*_data->load_cols.get()), *temp_ptr, ref_temp);
    // }
}

/**
 * Adds a support constraint to a named node region using the active support
 * collector and an optional coordinate system.
 */
void Model::add_support(const std::string& nset,
                        const StaticVector<6> constraint,
                        const std::string& orientation) {
    logging::error(_data->supp_cols.has_any(),
        "No support collectors have been defined");
    logging::error(_data->supp_cols.get() != nullptr,
        "No support collectors is currently active");

    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
    }

    bc::SupportCollector::Ptr supp_col = _data->supp_cols.get();
    supp_col->add_supp(
        _data->node_sets.get(nset),
        constraint,
        _data->coordinate_systems.get(orientation)
    );
}

/**
 * Adds a support constraint to one global node id using a temporary node region.
 */
void Model::add_support(ID id,
                        const StaticVector<6> constraint,
                        const std::string& orientation) {
    logging::error(_data->supp_cols.has_any(),
        "No support collectors have been defined");
    logging::error(_data->supp_cols.get() != nullptr,
        "No support collectors is currently active");

    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
    }

    NodeRegion::Ptr region = std::make_shared<NodeRegion>("INTERNAL");
    region->add(id);

    bc::SupportCollector::Ptr supp_col = _data->supp_cols.get();
    supp_col->add_supp(region, constraint, _data->coordinate_systems.get(orientation));
}

/**
 * Creates a solid section assignment for a named element set and material.
 * An optional coordinate system defines the material orientation.
 */
void Model::solid_section(const std::string& set,
                          const std::string& material,
                          const std::string& orientation) {
    logging::error(_data->elem_sets.has(set),
        "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material),
        "Material ", material, " is not a defined material");

    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
    }

    SolidSection::Ptr section = std::make_shared<SolidSection>();
    section->material_    = _data->materials.get(material);
    section->region_      = _data->elem_sets.get(set);
    section->orientation_ = orientation.empty() ? nullptr : _data->coordinate_systems.get(orientation);
    this->_data->sections.push_back(section);
}

/**
 * Creates a beam section assignment from an element set, material, profile and
 * beam-orientation direction.
 */
void Model::beam_section(const std::string& set,
                         const std::string& material,
                         const std::string& profile,
                         Vec3 orientation) {
    logging::error(_data->elem_sets.has(set),
        "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material),
        "Material ", material, " is not a defined material");
    logging::error(_data->profiles.has(profile),
        "Profile ", profile, " is not a defined profile");

    BeamSection::Ptr section = std::make_shared<BeamSection>();
    section->material_  = _data->materials.get(material);
    section->region_    = _data->elem_sets.get(set);
    section->profile_   = _data->profiles.get(profile);
    section->direction_ = orientation;
    this->_data->sections.push_back(section);
}

/**
 * Creates a truss section assignment with the prescribed cross-sectional area.
 */
void Model::truss_section(const std::string& set,
                          const std::string& material,
                          Precision area) {
    logging::error(_data->elem_sets.has(set),
        "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material),
        "Material ", material, " is not a defined material");

    TrussSection::Ptr section = std::make_shared<TrussSection>(
        _data->materials.get(material),
        _data->elem_sets.get(set),
        area
    );
    this->_data->sections.push_back(section);
}

/**
 * Creates an integrated shell section assignment.
 *
 * The section references one material, physical thickness, optional material
 * coordinate system and the selected local coordinate-system axis.
 */
void Model::shell_section(const std::string& set,
                          const std::string& material,
                          Precision thickness,
                          const std::string& orientation,
                          Index csys_axis) {
    logging::error(_data->elem_sets.has(set),
        "Element set ", set, " is not a defined element set");
    logging::error(_data->materials.has(material),
        "Material ", material, " is not a defined material");
    logging::error(csys_axis >= 0 && csys_axis < 3,
        "SHELLSECTION: CSYSAXIS must be 1, 2 or 3");

    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
    }

    IntegratedShellSection::Ptr section = std::make_shared<IntegratedShellSection>(
        _data->materials.get(material),
        _data->elem_sets.get(set),
        thickness,
        orientation.empty() ? nullptr : _data->coordinate_systems.get(orientation),
        csys_axis
    );
    this->_data->sections.push_back(section);
}

/**
 * Creates a shell section from prescribed ABD and transverse-shear matrices.
 *
 * The material is optional because the supplied constitutive matrices may define
 * the shell response directly. Region, thickness, optional orientation and local
 * coordinate-system axis are stored with the section.
 */
void Model::shell_section_abd(const std::string& set,
                              const std::string& material,
                              Precision thickness,
                              const Mat6& abd,
                              const Mat2& shear,
                              const std::string& orientation,
                              Index csys_axis) {
    logging::error(_data->elem_sets.has(set),
        "Element set ", set, " is not a defined element set");
    logging::error(csys_axis >= 0 && csys_axis < 3,
        "SHELLSECTION: CSYSAXIS must be 1, 2 or 3");

    if (!material.empty()) {
        logging::error(_data->materials.has(material),
            "Material ", material, " is not a defined material");
    }
    if (!orientation.empty()) {
        logging::error(_data->coordinate_systems.has(orientation),
            "Coordinate system ", orientation, " does not exist");
    }

    ABDShellSection::Ptr section = std::make_shared<ABDShellSection>(
        material.empty() ? nullptr : _data->materials.get(material),
        _data->elem_sets.get(set),
        thickness,
        abd,
        shear,
        orientation.empty() ? nullptr : _data->coordinate_systems.get(orientation),
        csys_axis
    );
    this->_data->sections.push_back(section);
}

/**
 * Adds a point-mass feature to a named node region.
 *
 * Translational mass, rotary inertia and translational/rotational spring
 * constants are stored on the feature and later contribute through model-level
 * mass and stiffness assembly.
 */
void Model::add_point_mass_feature(const std::string& nset,
                                   Precision mass,
                                   Vec3 rotary_inertia,
                                   Vec3 spring,
                                   Vec3 rotary_spring) {
    logging::error(_data->node_sets.has(nset),
        "Node set ", nset, " is not a defined node set");

    auto feature = std::make_shared<feature::PointMass>();
    feature->region_                  = _data->node_sets.get(nset);
    feature->mass_                    = mass;
    feature->rotary_inertia_          = rotary_inertia;
    feature->spring_constants_        = spring;
    feature->rotary_spring_constants_ = rotary_spring;
    this->_data->features.push_back(std::move(feature));
}

/**
 * Writes a compact model summary and emits detailed information for materials,
 * sections, profiles and named element/node sets through their existing `info()`
 * methods.
 *
 * @return The supplied output stream after the scalar capacity summary.
 */
std::ostream& operator<<(std::ostream& ostream, const model::Model& model) {
    ostream << "max nodes = " << model._data->max_nodes << '\n';
    ostream << "max elements = " << model._data->max_elems << '\n';
    ostream << "max surfaces = " << model._data->max_surfaces << "\n";

    logging::info(true, "Materials");
    logging::up();
    for (const auto& material : model._data->materials) {
        material.second->info();
    }
    logging::down();

    logging::info(true, "Sections");
    logging::up();
    for (const auto& section : model._data->sections) {
        section->info();
    }
    logging::down();

    logging::info(true, "Profiles");
    logging::up();
    for (const auto& profile : model._data->profiles) {
        profile.second->info();
    }
    logging::down();

    logging::info(true, "Element sets");
    logging::up();
    for (const auto& elem_set : model._data->elem_sets) {
        elem_set.second->info();
    }
    logging::down();

    logging::info(true, "Node sets");
    logging::up();
    for (const auto& node_set : model._data->node_sets) {
        node_set.second->info();
    }
    logging::down();

    return ostream;

    // // print materials
    // ostream << "Materials:\n";
    // for (const auto& material : model._data->materials) {
    //     ostream << *material.second;
    // }
    // return ostream;
}

} // namespace model
} // namespace fem
