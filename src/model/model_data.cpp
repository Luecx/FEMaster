/**
 * @file model_data.cpp
 * @brief Provides the translation unit companion for `ModelData`.
 *
 * Currently this file only includes the header to satisfy build systems that
 * expect a `.cpp` next to the declaration.
 *
 * @see src/model/model_data.h
 */

#include "model_data.h"
#include "element/element.h"

namespace fem {
namespace model {

Index ModelData::field_rows(FieldDomain domain) const {
    switch (domain) {
        case FieldDomain::UNKNOWN:
            logging::error(false, "ModelData: cannot allocate UNKNOWN fields");
            return 0;
        case FieldDomain::NODE:
            return static_cast<Index>(max_nodes);
        case FieldDomain::ELEMENT:
            return static_cast<Index>(max_elems);
        case FieldDomain::ELEMENT_NODAL:
            logging::error(element_nodal_offsets != nullptr,
                           "ModelData: element nodal offset field is not initialized");
            logging::error((*element_nodal_offsets)(static_cast<Index>(max_elems)) > 0,
                           "ModelData: no element nodes are available, cannot allocate ELEMENT_NODAL fields");
            return static_cast<Index>((*element_nodal_offsets)(static_cast<Index>(max_elems)));
        case FieldDomain::ELEMENT_IP:
            logging::error(element_ip_offsets != nullptr,
                           "ModelData: element IP offset field is not initialized");
            logging::error((*element_ip_offsets)(static_cast<Index>(max_elems)) > 0,
                           "ModelData: no integration points are available, cannot allocate ELEMENT_IP fields");
            return static_cast<Index>((*element_ip_offsets)(static_cast<Index>(max_elems)));
    }
    logging::error(false, "ModelData: unknown field domain");
    return 0;
}

void ModelData::register_element_offsets(ElementInterface& element) {
    logging::error(element_nodal_offsets != nullptr,
                   "ModelData: element nodal offset field is not initialized");
    logging::error(element_ip_offsets != nullptr,
                   "ModelData: element IP offset field is not initialized");
    logging::error(element.elem_id >= 0 && element.elem_id < max_elems,
                   "ModelData: element id ", element.elem_id, " out of bounds for offset fields");

    const Index first_affected = static_cast<Index>(element.elem_id + 1);
    const Index sentinel = static_cast<Index>(max_elems);
    const Precision node_delta = static_cast<Precision>(element.n_nodes());
    const Precision ip_delta = static_cast<Precision>(element.n_integration_points());

    for (Index row = first_affected; row <= sentinel; ++row) {
        (*element_nodal_offsets)(row) += node_delta;
        (*element_ip_offsets)(row) += ip_delta;
    }

    max_integration_points = static_cast<ID>((*element_ip_offsets)(sentinel));
}

bool ModelData::has_field(const std::string& name) const {
    return fields.find(name) != fields.end();
}

Field::Ptr ModelData::get_field(const std::string& name) const {
    auto it = fields.find(name);
    if (it == fields.end()) {
        return nullptr;
    }
    return it->second;
}

Field::Ptr ModelData::create_field(const std::string& name, FieldDomain domain, Index components, bool fill_nan, bool reg) {
    logging::error(!name.empty(), "Field name cannot be empty");

    if (!reg) {
        const Index rows = field_rows(domain);
        auto field = std::make_shared<Field>(name, domain, rows, components);
        if (fill_nan) {
            field->fill_nan();
        }
        return field;
    }

    auto it = fields.find(name);
    if (it != fields.end()) {
        Field::Ptr field = it->second;
        logging::error(field->domain == domain, "Field '", name, "': domain mismatch");
        logging::error(field->components == components, "Field '", name, "': components mismatch");
        return field;
    }

    const Index rows = field_rows(domain);
    auto field = std::make_shared<Field>(name, domain, rows, components);
    if (fill_nan) {
        field->fill_nan();
    }
    fields.emplace(name, field);
    return field;
}

Field ModelData::create_field_(const std::string& name, FieldDomain domain, Index components, bool fill_nan) {
    logging::error(!name.empty(), "Field name cannot be empty");

    const Index rows = field_rows(domain);
    auto field = Field(name, domain, rows, components);
    if (fill_nan) {
        field.fill_nan();
    }
    return field;
}
} // namespace model
} // namespace fem
