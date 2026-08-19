/**
 * @file model_data.h
 * @brief Declares semantic and compiled FEM model storage.
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
#include "../constraints/types/coupling.h"
#include "../constraints/types/contact.h"
#include "../constraints/types/equation.h"
#include "../constraints/types/rbm.h"
#include "../constraints/types/tie.h"
#include "../cos/coordinate_system.h"
#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../data/dict.h"
#include "../data/field.h"
#include "../data/region.h"
#include "../data/sets.h"
#include "../feature/feature.h"
#include "../material/material.h"
#include "../section/profile.h"
#include "../section/section.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace fem::model {

struct Part;
struct Instance;

struct ModelData {

    // non-flat data. Each part has some sets which should, once compiled, also be present in the sets below
    // instances reference a part. the sets will be copied for each instance.
    Dict<Part>     parts;
    Dict<Instance> instances;
    bool           compiled = false;

    // flattened elements, surfaces and lines. mapping to instance elements is
    // possible using the map inside each instance.
    std::vector<ElementPtr> elements;
    std::vector<SurfacePtr> surfaces;
    std::vector<LinePtr>    lines;

    // flattened sections, profiles and other features.
    std::vector<Section::Ptr>  sections;
    Dict<Profile>              profiles;
    std::vector<feature::Feature::Ptr> features;

    // materials, coordinate systems and amplitudes defined across the assembly. (there are no part-level csys)
    Dict<material::Material>    materials;
    Dict<cos::CoordinateSystem> coordinate_systems;
    Dict<bc::Amplitude>         amplitudes;

    // global fields specified for the full assembly.
    std::unordered_map<std::string, Field::Ptr> fields;

    // specific fields used throughout the solver. may be bound to a field from the fields map above.
    Field::Ptr positions                   = nullptr;
    Field::Ptr positions_reference         = nullptr;
    Field::Ptr element_stiffness_scale     = nullptr;
    Field::Ptr material_orientation        = nullptr;
    Field::Ptr material_state              = nullptr;
    Field::Ptr shell_element_nodal_normals = nullptr;
    Field::Ptr element_nodal_offsets       = nullptr;
    Field::Ptr element_ip_offsets          = nullptr;
    Field::Ptr element_mp_offsets          = nullptr;

    // sets for nodes, elements, surfaces and lines. should be just the sets from each instance / part.
    Sets<NodeRegion>    node_sets   {SET_NODE_ALL};
    Sets<ElementRegion> elem_sets   {SET_ELEM_ALL};
    Sets<SurfaceRegion> surface_sets{SET_SURF_ALL};
    Sets<LineRegion>    line_sets   {SET_LINE_ALL};

    // constraints defined for the assembly.
    std::vector<constraint::Connector> connectors;
    std::vector<constraint::Coupling>  couplings;
    std::vector<constraint::Tie>       ties;
    std::vector<constraint::Contact>   contacts;
    std::vector<constraint::Rbm>       rbms;
    std::vector<constraint::Equation>  equations;

    // all support collectors and load collectors defined for all the load-cases.
    Sets<bc::SupportCollector> supp_cols;
    Sets<bc::LoadCollector>    load_cols;

    // empty constructor
    ModelData() = default;

    // specifying how many total rows there are for a specific domain, e.g. Nodal, Element Nodal, and so on.
    Index field_rows(FieldDomain domain) const;
    // initialise the fields for the offsets for each field type for each element.
    void initialize_element_enumeration();

    // field creation and management
    bool has_field(const std::string& name) const;
    Field::Ptr get_field(const std::string& name) const;
    Field::Ptr create_field(const std::string& name, FieldDomain domain, Index components, bool fill_nan = true, bool reg = true);
    Field create_field_(const std::string& name, FieldDomain domain, Index components, bool fill_nan = true);

    // conversion of fields from one domain to another.
    Field element_nodal_to_nodal(const Field& element_nodal,
                                 const Field& element_weights,
                                 const std::string& name) const;
};

} // namespace fem::model
