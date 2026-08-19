/**
 * @file element.h
 * @brief Declares the base interface implemented by all finite elements.
 *
 * Part elements are persistent semantic definitions. During `Model::compile()`
 * each instance asks the concrete element for an independent polymorphic copy,
 * then rewires only assembly-owned ids and connectivity. The virtual `copy()`
 * function is the complete cloning contract; no external type registry or
 * stored function pointer is required.
 *
 * @see src/model/element/element_structural.h
 */

#pragma once

#include "../../core/types_cls.h"
#include "../../section/section.h"
#include "../geometry/surface/surface.h"
#include "../model_data.h"

#include <string>

namespace fem::model {

struct ElementInterface {
    ID elem_id           = 0;
    ID elem_nodal_offset = 0;
    ID elem_ip_offset    = 0;
    ID elem_mp_offset    = 0;

    Section::Ptr _section    = nullptr;
    ModelData*   _model_data = nullptr;

    explicit ElementInterface(ID elem_id_in)
        : elem_id(elem_id_in) {}

    virtual ~ElementInterface() = default;

    virtual ElDofs dofs() const       = 0;
    virtual Dim    dimensions() const = 0;
    virtual Dim    n_nodes() const    = 0;
    virtual Dim    num_ip() const     = 0;
    virtual const ID* nodes() const   = 0;

    /**
     * @brief Returns an independent copy preserving the concrete element type.
     *
     * The returned object represents the same persistent element definition as
     * the source. Runtime caches may intentionally be rebuilt instead of copied.
     * `Model::compile()` subsequently overwrites dense ids, offsets, section
     * assignment, ModelData binding and connectivity for the target Instance.
     */
    virtual ElementPtr copy() const = 0;

    virtual ID* nodes() {
        return const_cast<ID*>(static_cast<const ElementInterface*>(this)->nodes());
    }

    virtual Index num_mp_per_ip() const { return 1; }

    virtual SurfacePtr surface(ID) { return nullptr; }
    virtual LinePtr line(ID) { return nullptr; }

    virtual std::string type_name() const { return std::string{}; }

    ID* begin() { return nodes(); }
    ID* end() { return nodes() + n_nodes(); }
    const ID* begin() const { return nodes(); }
    const ID* end() const { return nodes() + n_nodes(); }

    Index ip_index(Index local_ip) const {
        return static_cast<Index>(elem_ip_offset) + local_ip;
    }

    Index mp_index(Index local_ip, Index local_mp = 0) const {
        return static_cast<Index>(elem_mp_offset)
             + local_ip * num_mp_per_ip()
             + local_mp;
    }

    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }

    template<typename T>
    const T* as() const {
        return dynamic_cast<const T*>(this);
    }

    void set_section(Section::Ptr section) {
        _section = std::move(section);
    }

    material::MaterialPtr material() {
        logging::error(_section != nullptr,
            "no section assigned to element ", elem_id);
        logging::error(_section->material_ != nullptr,
            "no material assigned to element ", elem_id);
        return _section->material_;
    }

    Vec3 node_position(ID local_id) const {
        logging::error(_model_data != nullptr,
            "no model data assigned to element ", elem_id);
        logging::error(_model_data->positions != nullptr,
            "positions field not set in model data");
        logging::error(local_id >= 0 && local_id < static_cast<ID>(n_nodes()),
            "local node id ", local_id, " out of range for element ", elem_id);
        return _model_data->positions->row_vec3(static_cast<Index>(nodes()[static_cast<Index>(local_id)]));
    }
};

} // namespace fem::model
