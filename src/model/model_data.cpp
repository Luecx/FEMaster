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

namespace fem {
namespace model {

Index ModelData::field_rows(FieldDomain domain) const {
    switch (domain) {
        case FieldDomain::NODE:
            return static_cast<Index>(max_nodes);
        case FieldDomain::ELEMENT:
            return static_cast<Index>(max_elems);
        case FieldDomain::IP:
            logging::error(max_integration_points > 0,
                           "ModelData: max_integration_points is 0, cannot allocate IP fields");
            return static_cast<Index>(max_integration_points);
    }
    logging::error(false, "ModelData: unknown field domain");
    return 0;
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
