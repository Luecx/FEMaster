/**
 * @file model.ipp
 * @brief Implements compact model construction operations.
 *
 * Topology construction writes directly into the active semantic Part. Once
 * `compile()` has flattened Parts through Instances, these operations are frozen
 * because changing local topology would invalidate dense assembly identifiers.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "element/element.h"
#include "element/element_structural.h"
#include "geometry/surface/surface.h"

inline void Model::set_node(ID id, Precision x, Precision y, Precision z) {
    logging::error(_data != nullptr && !_data->compiled,
        "Model: nodes cannot be added after compile()");

    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");
    logging::error(active->nodes.find(id) == active->nodes.end(),
        "Model: node ", id, " is already defined in part ", active->name);

    active->nodes.emplace(id, Vec3{x, y, z});
    active->node_sets.add(id);
}

template<typename T, typename... Args>
inline void Model::set_element(ID id, Args&&... args) {
    logging::error(_data != nullptr && !_data->compiled,
        "Model: elements cannot be added after compile()");

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

    active->elements.emplace(id, std::move(element));
    active->elem_sets.add(id);
}

inline void Model::set_surface(ID id, ID element_id, ID surface_id) {
    logging::error(_data != nullptr && !_data->compiled,
        "Model: surfaces cannot be added after compile()");

    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");

    const auto element_it = active->elements.find(element_id);
    logging::error(element_it != active->elements.end() && element_it->second != nullptr,
        "Model: element ", element_id, " is not defined in part ", active->name);

    auto surface = element_it->second->surface(surface_id);
    auto line    = element_it->second->line(surface_id);

    logging::error(surface != nullptr || line != nullptr,
        "Model: boundary ", surface_id, " of element ", element_id,
        " in part ", active->name, " provides neither a surface nor a line");

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

inline void Model::set_surface(const std::string& elset, ID surface_id) {
    logging::error(_data != nullptr && !_data->compiled,
        "Model: surfaces cannot be added after compile()");

    const auto active = _data->parts.get();
    logging::error(active != nullptr,
        "Model: no active part is available");
    logging::error(active->elem_sets.has(elset),
        "Model: element set ", elset, " is not defined in part ", active->name);
    logging::error(active->elem_sets.get(elset) && active->elem_sets.get(elset)->size() > 0,
        "Model: element set ", elset, " is empty in part ", active->name);

    for (const ID element_id : *active->elem_sets.get(elset)) {
        set_surface(-1, element_id, surface_id);
    }
}
