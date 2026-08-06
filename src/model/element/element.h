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
    const ID elem_id           = 0; ///< Unique element identifier.
    ID       elem_nodal_offset = 0; ///< Offset into element-nodal result fields.
    ID       elem_ip_offset    = 0; ///< Offset into integration-point result fields.
    ID       elem_mp_offset    = 0; ///< Offset into material-point result fields.

    Section::Ptr _section    = nullptr; ///< Associated section definition.
    ModelData*   _model_data = nullptr; ///< Non-owning back-reference to the owning model data.

    explicit ElementInterface(ID elem_id_in)
        : elem_id(elem_id_in) {}

    virtual ~ElementInterface() = default;

    virtual ElDofs    dofs() const                 = 0;
    virtual Dim       dimensions() const           = 0;
    virtual Dim       n_nodes() const              = 0;
    virtual Dim       num_ip() const               = 0;
    virtual const ID* nodes() const                = 0;

    virtual Index num_mp_per_ip() const { return 1; }

    virtual SurfacePtr surface(ID) { return nullptr; }
    virtual LinePtr    line   (ID) { return nullptr; }

    /// Short type tag (e.g., "C3D8", "S4", "B33"). Override in derived types.
    virtual std::string type_name() const { return std::string{}; }

    /// Iterator access over nodal identifiers.
    const ID* begin() const { return nodes(); }
    const ID* end  () const { return nodes() + n_nodes(); }

    Index ip_index(Index local_ip) const {
        return static_cast<Index>(elem_ip_offset) + local_ip;
    }

    Index mp_index(Index local_ip, Index local_mp = 0) const {
        return static_cast<Index>(elem_mp_offset)
             + local_ip * num_mp_per_ip()
             + local_mp;
    }

    /// Casts the element to a specific derived type, returning `nullptr` on failure.
    template<typename T> T* as() {
        return dynamic_cast<T*>(this);
    }

    /// Assigns the section used by the element.
    void set_section(Section::Ptr section) {
        _section = std::move(section);
    }

    /// Returns the material referenced by the assigned section.
    material::MaterialPtr material() {
        logging::error(_section != nullptr, "no section assigned to element ", elem_id);
        logging::error(_section->material_ != nullptr, "no material assigned to element ", elem_id);
        return _section->material_;
    }

    /// Returns the global position of an element-local node.
    Vec3 node_position(ID local_id) const {
        logging::error(_model_data            != nullptr, "no model data assigned to element ", elem_id);
        logging::error(_model_data->positions != nullptr, "positions field not set in model data");
        logging::error(local_id >= 0 && local_id < static_cast<ID>(n_nodes()),
                       "local node id ", local_id, " out of range for element ", elem_id);
        return _model_data->positions->row_vec3(static_cast<Index>(nodes()[static_cast<Index>(local_id)]));
    }
};
} // namespace model
} // namespace fem
