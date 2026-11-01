/**
 * @file model_data.cpp
 * @brief Implements field sizing and element-local enumeration.
 *
 * The implementation derives storage dimensions from compiled assembly data,
 * constructs prefix offsets for variable-size element domains, manages named and
 * temporary fields.
 *
 * Element formulations own the number and ordering of their local nodes,
 * integration points and material points. `ModelData` concatenates those local
 * ranges into model-wide fields and records the start offset on every element.
 *
 * @see ModelData
 * @see Field
 * @see ElementInterface
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "model_data.h"
#include "element/element.h"

namespace fem::model {

/**
 * Determines the number of rows required by one model field domain.
 *
 * Node and element domains derive directly from dense assembly storage.
 * Element-nodal, integration-point and material-point domains use the terminal
 * value of their compiled prefix-offset field, which equals the sum of the
 * corresponding local counts over all elements.
 *
 * Variable-size element domains therefore require
 * `initialize_element_enumeration()` to have completed. The unknown domain is
 * not allocatable and is rejected explicitly.
 *
 * @param domain Physical entity associated with each requested field row.
 * @return Total number of rows required for the domain.
 */
Index ModelData::field_rows(FieldDomain domain) const {
    // Resolve fixed-size and concatenated domains from their authoritative data
    const Index element_count = static_cast<Index>(elements.size());

    switch (domain) {
        case FieldDomain::UNKNOWN:
            logging::error(false,
                "ModelData: cannot allocate UNKNOWN fields");
            return 0;
        case FieldDomain::NODE:
            logging::error(positions != nullptr,
                "ModelData: POSITION field is not initialized");
            return positions->rows;
        case FieldDomain::ELEMENT:
            return element_count;
        case FieldDomain::ELEMENT_NODAL:
            logging::error(element_nodal_offsets != nullptr,
                "ModelData: element nodal offsets are not initialized");
            return static_cast<Index>((*element_nodal_offsets)(element_count));
        case FieldDomain::ELEMENT_IP:
            logging::error(element_ip_offsets != nullptr,
                "ModelData: element IP offsets are not initialized");
            return static_cast<Index>((*element_ip_offsets)(element_count));
        case FieldDomain::ELEMENT_MP:
            logging::error(element_mp_offsets != nullptr,
                "ModelData: element MP offsets are not initialized");
            return static_cast<Index>((*element_mp_offsets)(element_count));
    }

    logging::error(false,
        "ModelData: unknown field domain");
    return 0;
}

/**
 * Builds model-wide prefix enumerations for all element-local storage domains.
 *
 * Each offset field contains `element_count + 1` scalar entries. Entry `e`
 * stores the first row owned by element `e`, and the following entry closes its
 * half-open range. Null element slots retain an empty range. The final entry is
 * consequently the complete row count of the corresponding field domain.
 *
 * Every represented element receives its three starting offsets for direct
 * access during assembly and result recovery. Element-nodal ranges contain one
 * row per connected node, integration-point ranges contain `num_ip()` rows and
 * material-point ranges contain `num_ip() * num_mp_per_ip()` rows.
 *
 * Material-state fields are allocated lazily by Model::step_begin(), once
 * section assignments determine the necessary constitutive state width.
 * Enumeration is a one-time operation because changing offsets would invalidate
 * all dependent fields and element references.
 */
void ModelData::initialize_element_enumeration() {
    // Reject repeated initialization of offsets referenced by compiled elements
    logging::error(element_nodal_offsets == nullptr
                && element_ip_offsets    == nullptr
                && element_mp_offsets    == nullptr,
        "ModelData: element enumeration has already been initialized");

    // Allocate one prefix entry per element and a terminal total entry
    const Index element_count = static_cast<Index>(elements.size());
    const Index offset_rows   = element_count + 1;

    element_nodal_offsets = std::make_shared<Field>(
        "ELEMENT_NODAL_OFFSETS", FieldDomain::ELEMENT, offset_rows, 1);
    element_ip_offsets    = std::make_shared<Field>(
        "ELEMENT_IP_OFFSETS", FieldDomain::ELEMENT, offset_rows, 1);
    element_mp_offsets    = std::make_shared<Field>(
        "ELEMENT_MP_OFFSETS", FieldDomain::ELEMENT, offset_rows, 1);

    Index nodal_offset = 0;
    Index ip_offset    = 0;
    Index mp_offset    = 0;

    // Record every element range and publish its starting offsets in the element
    for (Index row = 0; row < element_count; ++row) {
        (*element_nodal_offsets)(row) = static_cast<Precision>(nodal_offset);
        (*element_ip_offsets)(row)    = static_cast<Precision>(ip_offset);
        (*element_mp_offsets)(row)    = static_cast<Precision>(mp_offset);

        const auto& element = elements[static_cast<std::size_t>(row)];
        if (!element) {
            continue;
        }

        element->elem_nodal_offset = static_cast<ID>(nodal_offset);
        element->elem_ip_offset    = static_cast<ID>(ip_offset);
        element->elem_mp_offset    = static_cast<ID>(mp_offset);

        nodal_offset += static_cast<Index>(element->n_nodes());
        ip_offset    += static_cast<Index>(element->num_ip());
        mp_offset    += static_cast<Index>(element->num_ip()) * element->num_mp_per_ip();
    }

    // Close the final element ranges and expose the total domain sizes
    (*element_nodal_offsets)(element_count) = static_cast<Precision>(nodal_offset);
    (*element_ip_offsets)(element_count)    = static_cast<Precision>(ip_offset);
    (*element_mp_offsets)(element_count)    = static_cast<Precision>(mp_offset);

    // Material history is initialized after sections are assigned, immediately
    // before the first element evaluation in Model::step_begin().
}

bool ModelData::has_field(const std::string& name) const {
    return fields.find(name) != fields.end();
}

Field::Ptr ModelData::get_field(const std::string& name) const {
    const auto it = fields.find(name);
    return it == fields.end() ? nullptr : it->second;
}

/**
 * Creates or retrieves a dense field with model-derived row dimensions.
 *
 * Registered creation first searches the global field dictionary. A matching
 * field is returned unchanged after its domain and component count have been
 * validated; `fill_nan` is not reapplied to existing storage. Otherwise a new
 * field is allocated using `field_rows(domain)`, optionally initialized with
 * NaN and inserted into the registry when requested.
 *
 * Unregistered fields remain independently owned by the returned shared pointer
 * and may use a diagnostic name already present in the registry.
 *
 * @param name Non-empty field name.
 * @param domain Physical entity associated with each field row.
 * @param components Number of scalar components stored per row.
 * @param fill_nan Initializes newly allocated storage with NaN when true.
 * @param reg Reuses and registers the field by name when true.
 * @return Shared pointer to the compatible existing or newly allocated field.
 */
Field::Ptr ModelData::create_field(const std::string& name,
                                   FieldDomain        domain,
                                   Index              components,
                                   bool               fill_nan,
                                   bool               reg) {
    // Validate the field identity before registry lookup or allocation
    logging::error(!name.empty(),
        "Field name cannot be empty");

    // Reuse registered storage only when its physical layout is compatible
    if (reg) {
        const auto it = fields.find(name);
        if (it != fields.end()) {
            const auto& field = it->second;
            logging::error(field->domain == domain,
                "Field '", name, "': domain mismatch");
            logging::error(field->components == components,
                "Field '", name, "': components mismatch");
            return field;
        }
    }

    // Allocate and initialize storage from the authoritative domain row count
    auto field = std::make_shared<Field>(name, domain, field_rows(domain), components);
    if (fill_nan) {
        field->fill_nan();
    }

    // Publish named persistent fields while leaving temporary fields independent
    if (reg) {
        fields.emplace(name, field);
    }
    return field;
}

/**
 * Creates an unregistered field value with model-derived row dimensions.
 *
 * This value-returning variant always allocates independent temporary storage
 * and never consults or modifies the named field registry.
 *
 * @param name Non-empty diagnostic field name.
 * @param domain Physical entity associated with each field row.
 * @param components Number of scalar components stored per row.
 * @param fill_nan Initializes the field with NaN when true.
 * @return Independently owned field value.
 */
Field ModelData::create_field_(const std::string& name,
                               FieldDomain        domain,
                               Index              components,
                               bool               fill_nan) const {
    // Validate and allocate temporary storage from the requested domain
    logging::error(!name.empty(),
        "Field name cannot be empty");

    Field field{name, domain, field_rows(domain), components};

    // Preserve zero initialization unless an undefined NaN field was requested
    if (fill_nan) {
        field.fill_nan();
    }
    return field;
}


} // namespace fem::model
