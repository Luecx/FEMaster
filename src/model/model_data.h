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
 * @struct ModelData
 * @brief Shared repository for model topology, fields, materials, and loads.
 */
struct ModelData {
    // Capacity information -----------------------------------------------------
    ID max_nodes;
    ID max_elems;
    ID max_surfaces;
    ID max_integration_points = 0;

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
    Field::Ptr positions = nullptr;
    Field::Ptr element_stiffness_scale = nullptr;
    Field::Ptr material_orientation = nullptr;

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
    ModelData(ID max_nodes, ID max_elems, ID max_surfaces, ID max_integration_points = 0)
        : max_nodes(max_nodes),
          max_elems(max_elems),
          max_surfaces(max_surfaces),
          max_integration_points(max_integration_points) {
        elements.resize(max_elems);
        surfaces.resize(max_surfaces);
    }

    // Field management --------------------------------------------------------

    /// Returns the row count associated with a field domain.
    Index field_rows(FieldDomain domain) const;

    /// Checks whether a named field exists.
    bool has_field(const std::string& name) const;

    /// Returns the field with the given name or nullptr if missing.
    Field::Ptr get_field(const std::string& name) const;

    /// Creates a field or returns the existing one after validation.
    Field::Ptr create_field (const std::string& name, FieldDomain domain, Index components, bool fill_nan = true, bool reg = true);
    Field      create_field_(const std::string& name, FieldDomain domain, Index components, bool fill_nan = true);
};

} // namespace model
} // namespace fem
