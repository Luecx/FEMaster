/**
 * @file model_data.h
 * @brief Declares the container that stores all FEM model input data.
 *
 * `ModelData` hosts elements, surfaces, materials, sections, fields, and the
 * various registries needed to assemble global matrices. It is shared across
 * higher-level model utilities and load-case builders.
 *
 * Element-local data is flattened into dense `ELEMENT_NODAL`, `ELEMENT_IP` and
 * `ELEMENT_MP` domains. Prefix-offset fields map each element to its contiguous
 * rows. Nonlinear constitutive updates expose committed and trial material-state
 * bindings explicitly through `material_state_old` and `material_state_new`.
 * `material_state` remains a compatibility alias to the active trial buffer for
 * elasticity-only element paths that have not yet migrated to the new contract.
 *
 * @see src/model/model_data.cpp
 * @see src/model/model.h
 *
 * @author Finn Eggers
 * @date 07.08.2026
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
#include "../constraints/types/tie.h"
#include "../constraints/types/rbm.h"
#include "../material/material.h"
#include "../cos/coordinate_system.h"
#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../data/dict.h"
#include "../data/field.h"
#include "../data/region.h"
#include "../data/sets.h"
#include "../section/profile.h"
#include "../section/section.h"
#include "../feature/feature.h"

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem {
namespace model {

/**
 * @brief Shared repository for model topology, fields and named definitions.
 *
 * `ModelData` owns the containers referenced by elements, sections, load cases
 * and writers. Entity arrays use stable global identifiers, while generic
 * fields carry their own semantic domain and component count.
 *
 * After all elements have been read, `initialize_element_enumeration()` assigns
 * contiguous element-nodal, integration-point and material-point spans. The
 * final sentinel entry of each offset field stores the total number of rows.
 * These offsets remain invariant for the lifetime of the assembled topology.
 *
 * Cached semantic field pointers provide fast access to positions, optional
 * scaling/orientation data and the active constitutive state bindings. A
 * nonlinear solve may rebind committed/trial material state, but field
 * dimensions must continue to match the enumerated `ELEMENT_MP` domain.
 */
struct ModelData {
    // Capacity information -----------------------------------------------------
    ID max_nodes;
    ID max_elems;
    ID max_surfaces;
    ID max_integration_points = 0;
    ID max_material_points    = 0;

    // Geometric entities -------------------------------------------------------
    std::vector<ElementPtr> elements;
    std::vector<SurfacePtr> surfaces;
    std::vector<LinePtr>    lines;

    // Sections and profiles ----------------------------------------------------
    std::vector<Section::Ptr> sections;
    Dict<Profile> profiles;

    // Non-element features ----------------------------------------------------
    std::vector<feature::Feature::Ptr> features;

    // Generic fields -----------------------------------------------------------
    std::unordered_map<std::string, Field::Ptr> fields;

    // Cached semantic fields ---------------------------------------------------
    Field::Ptr positions                   = nullptr;
    Field::Ptr positions_reference         = nullptr;
    Field::Ptr element_stiffness_scale     = nullptr;
    Field::Ptr material_orientation        = nullptr;
    Field::Ptr material_state_old          = nullptr;
    Field::Ptr material_state_new          = nullptr;
    Field::Ptr material_state              = nullptr; // compatibility alias to trial state
    Field::Ptr shell_element_nodal_normals = nullptr;
    Field::Ptr element_nodal_offsets       = nullptr;
    Field::Ptr element_ip_offsets          = nullptr;
    Field::Ptr element_mp_offsets          = nullptr;

    // Nonlinear kinematics -----------------------------------------------------
    // Existing nonlinear analyses remain geometrically nonlinear by default.
    // NLGEOM=OFF keeps constitutive nonlinearity but selects infinitesimal solid
    // kinematics and suppresses geometric stiffness.
    bool geometric_nonlinearity = true;

    // Region registries --------------------------------------------------------
    Sets<NodeRegion   > node_sets   {SET_NODE_ALL};
    Sets<ElementRegion> elem_sets   {SET_ELEM_ALL};
    Sets<SurfaceRegion> surface_sets{SET_SURF_ALL};
    Sets<LineRegion   > line_sets   {SET_LINE_ALL};

    // Named resources ----------------------------------------------------------
    Dict<material::Material> materials;
    Dict<cos::CoordinateSystem> coordinate_systems;
    Dict<bc::Amplitude> amplitudes;

    // Constraints --------------------------------------------------------------
    std::vector<constraint::Connector> connectors{};
    std::vector<constraint::Coupling> couplings{};
    std::vector<constraint::Tie> ties{};
    std::vector<constraint::Contact> contacts{};
    std::vector<constraint::Rbm> rbms{};
    std::vector<constraint::Equation> equations{};

    // Load and support collectors ---------------------------------------------
    Sets<bc::SupportCollector> supp_cols{};
    Sets<bc::LoadCollector> load_cols{};

    // Construction with fixed entity capacities. Element-local offset fields
    // remain uninitialized until the complete topology has been read.
    ModelData(ID max_nodes, ID max_elems, ID max_surfaces, ID max_integration_points = 0)
        : max_nodes(max_nodes),
          max_elems(max_elems),
          max_surfaces(max_surfaces),
          max_integration_points(max_integration_points),
          max_material_points(max_integration_points) {
        elements.resize(max_elems);
        surfaces.resize(max_surfaces);

        element_nodal_offsets = nullptr;
        element_ip_offsets    = nullptr;
        element_mp_offsets    = nullptr;
    }

    // Field sizing and element-local enumeration. Offset-dependent domains are
    // available only after initialize_element_enumeration() has completed.
    Index field_rows(FieldDomain domain);
    void initialize_element_enumeration();

    // Named field lookup and allocation. Registered creation reuses a compatible
    // field with the same name; unregistered creation returns temporary storage.
    bool has_field(const std::string& name) const;
    Field::Ptr get_field(const std::string& name) const;
    Field::Ptr create_field (const std::string& name, FieldDomain domain, Index components, bool fill_nan = true, bool reg = true);
    Field      create_field_(const std::string& name, FieldDomain domain, Index components, bool fill_nan = true);

    // Weighted projection from flattened element-nodal rows to shared nodes
    Field element_nodal_to_nodal(const Field& element_nodal,
                                 const Field& element_weights,
                                 const std::string& name) const;
};
} // namespace model
} // namespace fem
