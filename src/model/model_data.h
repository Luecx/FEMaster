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
 * The data root also owns shared material/profile definitions, compiled section
 * assignments, global field registration and prefix-sum enumerations for
 * element-nodal, integration-point and material-point storage. Public model
 * operations coordinate transitions and validation; `ModelData` provides the
 * common storage and domain-level utilities.
 *
 * @see Model
 * @see Part
 * @see Instance
 * @see Field
 *
 * @author Finn Eggers
 * @date 19.08.2026
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
 * `ModelData` is the internal ownership root behind `Model`. During model input,
 * Parts contain sparse local topology and Instances define rigid placements of
 * those reusable definitions. Compilation preserves this semantic data for name
 * and identifier resolution while additionally creating dense global arrays for
 * nodes, elements, surfaces, lines, regions and section assignments.
 *
 * Dense element identifiers index `elements` directly. Every represented
 * element stores the corresponding global nodal connectivity and receives
 * prefix offsets into the model-wide `ELEMENT_NODAL`, `ELEMENT_IP` and
 * `ELEMENT_MP` domains. Each offset field has one terminal entry so the span of
 * element `e` is `[offset[e], offset[e + 1])`.
 *
 * Named fields are owned through shared pointers in `fields`. Dedicated handles
 * expose fields with central solver meaning, such as current/reference nodal
 * positions, material history and shell reference normals. A dedicated handle
 * may also refer to a registered named field; the registry and handle therefore
 * share ownership rather than duplicate values.
 *
 * The `compiled` flag is monotonic. Once set, semantic topology must no longer be
 * changed because all dense identifiers, regions, offsets and dependent fields
 * rely on the compiled ordering.
 */
struct ModelData {

    // Semantic topology retained across compilation. Parts own sparse local
    // entities and regions; instances reference a Part and define its rigid
    // placement. `compiled` marks the one-way transition to dense assembly data.
    Dict<Part>     parts;
    Dict<Instance> instances;
    bool           compiled = false;

    // mapping global node back to local instance nodes
    std::vector<std::tuple<Instance::Ptr, Index>> node_mapping;

    // Dense assembly topology produced by Model::compile(). Vector indices are
    // global identifiers, while each Instance retains the corresponding map from
    // its part-local identifiers.
    std::vector<ElementPtr> elements;
    std::vector<SurfacePtr> surfaces;
    std::vector<LinePtr>    lines;

    // Assembly section assignments and shared non-topological definitions.
    // Features contribute matrices or loads outside regular element topology.
    std::vector<Section::Ptr>          sections;
    Dict<Profile>                      profiles;
    std::vector<feature::Feature::Ptr> features;

    // Named definitions shared across the complete assembly. Coordinate systems
    // and amplitudes are global resources rather than part-local topology.
    Dict<material::Material>    materials;
    Dict<cos::CoordinateSystem> coordinate_systems;
    Dict<bc::Amplitude>         amplitudes;

    // Named global fields indexed by their physical FieldDomain
    std::unordered_map<std::string, Field::Ptr> fields;

    // Solver-facing field handles. These may alias entries in `fields`; the
    // element offset fields are internal prefix sums with one terminal row.
    Field::Ptr positions                   = nullptr;
    Field::Ptr positions_reference         = nullptr;
    Field::Ptr element_stiffness_scale     = nullptr;
    Field::Ptr material_orientation        = nullptr;
    Field::Ptr material_state              = nullptr;
    Field::Ptr shell_element_nodal_normals = nullptr;
    Field::Ptr element_nodal_offsets       = nullptr;
    Field::Ptr element_ip_offsets          = nullptr;
    Field::Ptr element_mp_offsets          = nullptr;

    // Global assembly regions materialized from the compiled instances. Every
    // collection owns an aggregate *ALL region that receives all dense entities.
    Sets<NodeRegion>    node_sets   {SET_NODE_ALL};
    Sets<ElementRegion> elem_sets   {SET_ELEM_ALL};
    Sets<SurfaceRegion> surface_sets{SET_SURF_ALL};
    Sets<LineRegion>    line_sets   {SET_LINE_ALL};

    // Assembly constraints operating on dense global identifiers and regions
    std::vector<constraint::Connector> connectors;
    std::vector<constraint::Coupling>  couplings;
    std::vector<constraint::Tie>       ties;
    std::vector<constraint::Contact>   contacts;
    std::vector<constraint::Rbm>       rbms;
    std::vector<constraint::Equation>  equations;

    // Named support and load collectors shared by the model load cases
    Sets<bc::SupportCollector> supp_cols;
    Sets<bc::LoadCollector>    load_cols;

    // Construction
    ModelData() = default;

    // Domain sizing and element-local enumeration. Variable-size element domains
    // use terminal prefix offsets, whereas NODE and ELEMENT derive their row
    // counts directly from compiled coordinate and topology storage.
    Index field_rows(FieldDomain domain) const;
    void initialize_element_enumeration();

    // Named field lookup and allocation. Registered creation reuses a compatible
    // existing field; unregistered creation produces independent temporary
    // storage without changing the global field dictionary.
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

    // Weighted projection from element-local nodal rows onto unique global
    // nodes. Element weights control participation and the relative contribution
    // of every connected element to the final nodal average.
    Field element_nodal_to_nodal(
        const Field&       element_nodal,
        const Field&       element_weights,
        const std::string& name
    ) const;
};

} // namespace fem::model
