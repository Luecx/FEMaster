/**
 * @file model_data.h
 * @brief Declares the container that stores all FEM model input data.
 *
 * `ModelData` hosts elements, surfaces, materials, sections, fields, and the
 * various registries needed to assemble global matrices. It is shared across
 * higher-level model utilities and load-case builders.
 *
 * @see src/model/model_data.cpp
 * @see src/model/model.h
 */

#pragma once

#include "../bc/amplitude.h"
#include "../bc/load.h"
#include "../bc/load_collector.h"
#include "../bc/support.h"
#include "../bc/support_collector.h"
#include "../constraints/connector.h"
#include "../constraints/coupling.h"
#include "../constraints/equation.h"
#include "../constraints/tie.h"
#include "../material/material.h"
#include "../cos/coordinate_system.h"
#include "../core/types_eig.h"
#include "../data/dict.h"
#include "../data/elem_data_dict.h"
#include "../data/ip_data_dict.h"
#include "../data/node_data_dict.h"
#include "../data/region.h"
#include "../data/sets.h"
#include "../section/profile.h"
#include "../section/section.h"

#include <memory>
#include <vector>
#include <array>

namespace fem {
namespace model {

/**
 * @struct ModelData
 * @brief Shared repository for model topology, fields, materials, and loads.
 */
struct ModelData {
    // Capacity information -----------------------------------------------------
    ID max_nodes;
    ID max_elems;
    ID max_surfaces;

    // Geometric entities -------------------------------------------------------
    std::vector<ElementPtr> elements;
    std::vector<SurfacePtr> surfaces;
    // For 1D geometry (lines) used by line-based ties/couplings
    std::vector<std::array<ID, 2>> lines; ///< Currently supports 2-node lines (e.g., beam axes)

    // Sections and profiles ----------------------------------------------------
    std::vector<Section::Ptr> sections;
    Dict<Profile> profiles;

    // Core data tables ---------------------------------------------------------
    NodeDataDict node_data;
    ElemDataDict elem_data;

    // Field dictionaries -------------------------------------------------------
    NodeFieldDict node_fields;
    ElemFieldDict elem_fields;

    // Region registries --------------------------------------------------------
    Sets<NodeRegion> node_sets{SET_NODE_ALL};
    Sets<ElementRegion> elem_sets{SET_ELEM_ALL};
    Sets<SurfaceRegion> surface_sets{SET_SURF_ALL};
    Sets<LineRegion> line_sets{SET_LINE_ALL};

    // Named resources ----------------------------------------------------------
    Dict<material::Material> materials;
    Dict<cos::CoordinateSystem> coordinate_systems;
    Dict<bc::Amplitude> amplitudes;

    // Constraints --------------------------------------------------------------
    std::vector<constraint::Connector> connectors{};
    std::vector<constraint::Coupling> couplings{};
    std::vector<constraint::Tie> ties{};
    std::vector<constraint::Equation> equations{};

    // Load and support collectors ---------------------------------------------
    Sets<bc::SupportCollector> supp_cols{};
    Sets<bc::LoadCollector> load_cols{};

    /**
     * @brief Constructs the data repository with preallocated containers.
     */
    ModelData(ID max_nodes, ID max_elems, ID max_surfaces)
        : max_nodes(max_nodes),
          max_elems(max_elems),
          max_surfaces(max_surfaces),
          node_data(max_nodes),
          elem_data(max_elems) {
        elements.resize(max_elems);
        surfaces.resize(max_surfaces);
    }

    // Data management ---------------------------------------------------------

    /// Creates nodal data storage for the specified entry.
    void create_data(NodeDataEntries key, int entries) {
        node_data.create(key, entries);
    }

    /// Creates element data storage for the specified entry.
    void create_data(ElementDataEntries key, int entries) {
        elem_data.create(key, entries);
    }

    /// Returns mutable access to the nodal data block referenced by `key`.
    NodeData& get(NodeDataEntries key) {
        return node_data.get(key);
    }

    /// Returns mutable access to the element data block referenced by `key`.
    ElementData& get(ElementDataEntries key) {
        return elem_data.get(key);
    }

    /// Removes nodal data for the given entry.
    void remove(NodeDataEntries key) {
        node_data.remove(key);
    }

    /// Removes element data for the given entry.
    void remove(ElementDataEntries key) {
        elem_data.remove(key);
    }
};

} // namespace model
} // namespace fem
