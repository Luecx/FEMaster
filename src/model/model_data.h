/**
 * @file model_data.h
 * @brief Declares semantic and compiled FEM model storage.
 *
 * `ModelData` is the shared data root owned by `Model`. Before compilation it
 * retains reusable Parts and their rigid Instances. `Model::compile()` flattens
 * that semantic topology into dense assembly arrays, global regions and nodal
 * coordinate fields consumed by element formulations, constraints, load cases
 * and result writers.
 *
 * @see Model
 * @see Part
 * @see Instance
 * @see Field
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include "../bc/amplitude.h"
#include "../bc/load.h"
#include "../bc/load_collector.h"
#include "../bc/support.h"
#include "../bc/support_collector.h"
#include "../constraints/types/connector.h"
#include "../constraints/types/contact.h"
#include "../constraints/types/coupling.h"
#include "../constraints/types/equation.h"
#include "../constraints/types/rbm.h"
#include "../constraints/types/tie.h"
#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"
#include "../data/dict.h"
#include "../data/field.h"
#include "../data/region.h"
#include "../data/sets.h"
#include "../feature/feature.h"
#include "../material/material.h"
#include "../section/profile.h"
#include "../section/section.h"
#include "instance.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace fem::model {

struct Part;
struct Instance;

/**
 * @brief Shared semantic, assembly and solver-facing storage of a FEM model.
 *
 * Named fields are owned through shared pointers in `fields`. Dedicated handles
 * expose fields with central solver meaning. Constitutive history deliberately
 * uses two handles during nonlinear analysis: `material_state` is the immutable
 * source history for one physical increment and `material_state_target` is the
 * separate destination written by constitutive integration. Read-only result
 * recovery sees only the source state and therefore cannot advance history.
 */
struct ModelData {

    Dict<Part>     parts;
    Dict<Instance> instances;
    bool           compiled = false;

    // mapping global node back to local instance nodes
    std::vector<std::tuple<Instance::Ptr, Index>> node_mapping;

    std::vector<ElementPtr> elements;
    std::vector<SurfacePtr> surfaces;
    std::vector<LinePtr>    lines;

    std::vector<ElementPtr> point_elements;

    std::vector<Section::Ptr>          sections;
    Dict<Profile>                      profiles;
    std::vector<feature::Feature::Ptr> features;

    Dict<material::Material>    materials;
    Dict<cos::CoordinateSystem> coordinate_systems;
    Dict<bc::Amplitude>         amplitudes;

    std::unordered_map<std::string, Field::Ptr> fields;

    // Solver-facing field handles. `material_state` is always a source state;
    // stateful constitutive integration writes exclusively to
    // `material_state_target`. Both use identical ELEMENT_MP enumeration and
    // component layout when they are present.
    Field::Ptr positions                   = nullptr;
    Field::Ptr positions_reference         = nullptr;
    Field::Ptr element_stiffness_scale     = nullptr;
    Field::Ptr material_orientation        = nullptr;
    Field::Ptr material_state              = nullptr;
    Field::Ptr material_state_target       = nullptr;
    Field::Ptr shell_element_nodal_normals = nullptr;
    Field::Ptr element_nodal_offsets       = nullptr;
    Field::Ptr element_ip_offsets          = nullptr;
    Field::Ptr element_mp_offsets          = nullptr;

    Sets<NodeRegion>    node_sets   {SET_NODE_ALL};
    Sets<ElementRegion> elem_sets   {SET_ELEM_ALL};
    Sets<SurfaceRegion> surface_sets{SET_SURF_ALL};
    Sets<LineRegion>    line_sets   {SET_LINE_ALL};

    std::vector<constraint::Connector> connectors;
    std::vector<constraint::Coupling>  couplings;
    std::vector<constraint::Tie>       ties;
    std::vector<constraint::Contact>   contacts;
    std::vector<constraint::Rbm>       rbms;
    std::vector<constraint::Equation>  equations;

    Sets<bc::SupportCollector> supp_cols;
    Sets<bc::LoadCollector>    load_cols;

    ModelData() = default;

    Index field_rows(FieldDomain domain) const;
    void initialize_element_enumeration();

    bool has_field(const std::string& name) const;
    Field::Ptr get_field(const std::string& name) const;
    Field::Ptr create_field(
        const std::string& name,
        FieldDomain       domain,
        Index             components,
        bool              fill_nan = true,
        bool              reg = true
    );
    Field create_field_(
        const std::string& name,
        FieldDomain       domain,
        Index             components,
        bool              fill_nan = true
    );

    Field element_nodal_to_nodal(
        const Field&       element_nodal,
        const Field&       element_weights,
        const std::string& name
    ) const;
};

} // namespace fem::model