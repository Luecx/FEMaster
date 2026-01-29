/**
 * @file element.h
 * @brief Declares the base interface implemented by all finite elements.
 *
 * `ElementInterface` provides accessors shared across structural and non-
 * structural elements such as connectivity, sections, and model data.
 *
 * @see src/model/element/element_structural.h
 */

#pragma once

#include "../../core/types_cls.h"
#include <string>
#include "../../section/section.h"
#include "../geometry/surface/surface.h"
#include "../model_data.h"

namespace fem {
namespace model {

/**
 * @struct ElementInterface
 * @brief Minimal polymorphic interface for FEM elements.
 */
struct ElementInterface {
    const ID elem_id = 0; ///< Unique element identifier.

    Section::Ptr _section = nullptr; ///< Associated section definition.
    ModelDataPtr _model_data;        ///< Back-reference to the owning model data.

    explicit ElementInterface(ID elem_id_in)
        : elem_id(elem_id_in) {}

    virtual ~ElementInterface() = default;

    /// Assigns the section used by the element.
    void set_section(Section::Ptr section) { _section = std::move(section); }

    /// Returns the material referenced by the assigned section.
    material::MaterialPtr material() {
        logging::error(_section != nullptr, "no section assigned to element ", elem_id);
        logging::error(_section->material != nullptr, "no material assigned to element ", elem_id);
        return _section->material;
    }

    /// Returns the global position of a node from the positions field.
    Vec3 node_position(ID node_id) const {
        logging::error(_model_data != nullptr, "no model data assigned to element ", elem_id);
        logging::error(_model_data->positions != nullptr, "positions field not set in model data");
        const auto& positions = *_model_data->positions;
        const Index row = static_cast<Index>(node_id);
        return positions.row_vec3(row);
    }

    virtual ElDofs dofs() = 0;
    virtual Dim dimensions() = 0;
    virtual Dim n_nodes() = 0;
    virtual Dim n_integration_points() = 0;
    virtual ID* nodes() = 0;

    virtual SurfacePtr surface(ID surface_id) {return nullptr;};
    virtual LinePtr    line   (ID surface_id) {return nullptr;};

    /// Short type tag (e.g., "C3D8", "S4", "B33"). Override in derived types.
    virtual std::string type_name() const { return std::string{}; }

    /// Iterator access over nodal identifiers.
    ID* begin() { return nodes(); }
    ID* end() { return nodes() + n_nodes(); }

    /// Casts the element to a specific derived type, returning `nullptr` on failure.
    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }
};

} // namespace model
} // namespace fem
