/**
 * @file model.ipp
 * @brief Implements compact model construction operations.
 *
 * Topology construction writes directly into the active semantic Part. Once
 * `compile()` has flattened Parts through Instances, these operations are frozen
 * because changing local topology would invalidate dense assembly identifiers.
 *
 * Nodes retain part-local coordinates and sparse input identifiers. Element
 * construction stores local connectivity through the polymorphic element
 * interface, while boundary extraction asks the concrete element topology for
 * its corresponding surface or line. Every inserted entity is added to the
 * currently active local region when present and always to its aggregate *ALL
 * region.
 *
 * @see Model
 * @see Part
 * @see Model::compile
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "element/element.h"
#include "element/element_structural.h"
#include "geometry/surface/surface.h"

/**
 * Adds one node to the currently active semantic part.
 *
 * The supplied coordinates are stored in the local part coordinate system and
 * are transformed only when an Instance is compiled into the global assembly.
 * Node identifiers may be sparse but must be unique within the active part. The
 * new identifier is added to the active node set when present and always to the
 * part-local aggregate `NALL` region.
 *
 * @param id Part-local node identifier.
 * @param x Local x-coordinate.
 * @param y Local y-coordinate.
 * @param z Local z-coordinate.
 */
inline void Model::set_node(ID id, Precision x, Precision y, Precision z) {
    // Validate that semantic topology is still mutable
    logging::error(_data != nullptr && !_data->compiled,
        "Model: nodes cannot be added after compile()");

    // Resolve the active Part and enforce uniqueness in its local node namespace
    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");
    logging::error(active->nodes.find(id) == active->nodes.end(),
        "Model: node ", id, " is already defined in part ", active->name);

    // Store the local coordinates and register the identifier in active regions
    active->nodes.emplace(id, Vec3{x, y, z});
    active->node_sets.add(id);
}

/**
 * Constructs one concrete element in the currently active semantic part.
 *
 * Connectivity arguments are converted to part-local node identifiers and
 * passed to the element constructor in their supplied topology-specific order.
 * Connectivity existence is validated later by `Model::compile()`, when the
 * element is copied and rewired to dense assembly identifiers. The local
 * identifier is added to the active element set when present and always to the
 * aggregate `EALL` region.
 *
 * @tparam T Concrete element type derived from `ElementInterface`.
 * @tparam Args Node-identifier argument types accepted by the caller.
 * @param id Part-local element identifier.
 * @param args Ordered part-local nodal connectivity.
 */
template<typename T, typename... Args>
inline void Model::set_element(ID id, Args&&... args) {
    // Validate that semantic topology is still mutable
    logging::error(_data != nullptr && !_data->compiled,
        "Model: elements cannot be added after compile()");

    // Resolve the active Part and enforce uniqueness in its element namespace
    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");
    logging::error(active->elements.find(id) == active->elements.end(),
        "Model: element ", id, " is already defined in part ", active->name);

    // The concrete element owns its polymorphic copy() implementation. No
    // Model-side type registry is required; compile() can clone through the
    // ElementInterface pointer and then rewire the returned instance.
    auto element = std::make_shared<T>(id, std::array<ID, sizeof...(Args)>{
        static_cast<ID>(std::forward<Args>(args))...
    });

    // Transfer the element into local topology and register its identifier
    active->elements.emplace(id, std::move(element));
    active->elem_sets.add(id);
}

/**
 * Extracts one topological boundary from an element in the active part.
 *
 * The concrete element decides whether the requested boundary identifier
 * represents a finite-element surface, a line, or both. Extracted connectivity
 * remains part-local until compilation rewires it through each Instance. A
 * non-negative `id` is used directly; a negative value selects the lowest free
 * identifier independently in the surface and line namespaces.
 *
 * Every created boundary is inserted into the corresponding active region when
 * present and always into the aggregate `SFALL` or `LALL` region.
 *
 * @param id Requested part-local boundary identifier, or a negative value for
 *           automatic allocation.
 * @param element_id Part-local identifier of the source element.
 * @param surface_id Element-specific boundary identifier.
 */
inline void Model::set_surface(ID id, ID element_id, ID surface_id) {
    // Validate that semantic topology is still mutable
    logging::error(_data != nullptr && !_data->compiled,
        "Model: surfaces cannot be added after compile()");

    // Resolve the active Part and the element that owns the requested boundary
    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");

    const auto element_it = active->elements.find(element_id);
    logging::error(element_it != active->elements.end() && element_it->second != nullptr,
        "Model: element ", element_id, " is not defined in part ", active->name);

    // Ask the concrete element topology for its surface or line representation
    auto surface = element_it->second->surface(surface_id);
    auto line    = element_it->second->line(surface_id);

    logging::error(surface != nullptr || line != nullptr,
        "Model: boundary ", surface_id, " of element ", element_id,
        " in part ", active->name, " provides neither a surface nor a line");

    // Allocate and register a part-local surface when the element provides one
    if (surface) {
        ID local_id = id;
        if (local_id < 0) {
            local_id = 0;
            while (active->surfaces.find(local_id) != active->surfaces.end()) {
                ++local_id;
            }
        }

        logging::error(active->surfaces.find(local_id) == active->surfaces.end(),
            "Model: surface ", local_id, " is already defined in part ", active->name);

        active->surfaces.emplace(local_id, std::move(surface));
        active->surface_sets.add(local_id);
    }

    // Allocate and register a part-local line when the element provides one
    if (line) {
        ID local_id = id;
        if (local_id < 0) {
            local_id = 0;
            while (active->lines.find(local_id) != active->lines.end()) {
                ++local_id;
            }
        }

        logging::error(active->lines.find(local_id) == active->lines.end(),
            "Model: line ", local_id, " is already defined in part ", active->name);

        active->lines.emplace(local_id, std::move(line));
        active->line_sets.add(local_id);
    }
}

/**
 * Extracts the same element boundary for every element in a local region.
 *
 * The named element region must belong to the active part and contain at least
 * one element. Each boundary is delegated to the single-element overload with
 * automatic identifier allocation, selecting the lowest currently free
 * identifier in the applicable local surface or line namespace.
 *
 * @param elset Name of the part-local element region to expand.
 * @param surface_id Element-specific boundary identifier extracted from every
 *                   selected element.
 */
inline void Model::set_surface(const std::string& elset, ID surface_id) {
    // Validate that semantic topology is still mutable
    logging::error(_data != nullptr && !_data->compiled,
        "Model: surfaces cannot be added after compile()");

    // Resolve and validate the non-empty element region in the active Part
    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");
    logging::error(active->elem_sets.has(elset),
        "Model: element set ", elset, " is not defined in part ", active->name);
    logging::error(active->elem_sets.get(elset) && active->elem_sets.get(elset)->size() > 0,
        "Model: element set ", elset, " is empty in part ", active->name);

    // Extract one automatically numbered boundary from every selected element
    for (const ID element_id : *active->elem_sets.get(elset)) {
        set_surface(-1, element_id, surface_id);
    }
}
