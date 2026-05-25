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

Index ModelData::field_rows(FieldDomain domain) {
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

void ModelData::initialize_element_enumeration() {
    logging::error(element_nodal_offsets == nullptr && element_ip_offsets == nullptr,
                   "ModelData: element enumeration has already been initialized");

    element_nodal_offsets = std::make_shared<Field>(
        "ELEMENT_NODAL_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);
    element_ip_offsets = std::make_shared<Field>(
        "ELEMENT_IP_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);

    Index nodal_offset = 0;
    Index ip_offset = 0;

    for (Index row = 0; row < static_cast<Index>(max_elems); ++row) {
        (*element_nodal_offsets)(row) = static_cast<Precision>(nodal_offset);
        (*element_ip_offsets)(row) = static_cast<Precision>(ip_offset);

        const ElementPtr& element = elements[row];
        if (element != nullptr) {
            nodal_offset += static_cast<Index>(element->n_nodes());
            ip_offset += static_cast<Index>(element->n_integration_points());
        }
    }

    const Index sentinel = static_cast<Index>(max_elems);
    (*element_nodal_offsets)(sentinel) = static_cast<Precision>(nodal_offset);
    (*element_ip_offsets)(sentinel) = static_cast<Precision>(ip_offset);
    max_integration_points = static_cast<ID>(ip_offset);
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
