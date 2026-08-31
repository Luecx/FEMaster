/**
 * @file model_data.cpp
 * @brief Implements field sizing and element-local enumeration.
 *
 * The implementation derives storage dimensions from compiled assembly data,
 * constructs prefix offsets for variable-size element domains, manages named and
 * temporary fields and projects element-nodal results onto unique global nodes.
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
 * When the model contains material points, a zeroed single-component state field
 * is installed as the solver-facing default. Nonlinear state management may
 * replace it with a wider constitutive history field before material evaluation.
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

    // Provide separate zeroed input and output state rows for every enumerated
    // material point. Constitutive evaluations may never overwrite their input.
    if (mp_offset > 0) {
        material_state_old = std::make_shared<Field>(
            "MATERIAL_STATE_OLD", FieldDomain::ELEMENT_MP, mp_offset, 1);
        material_state_new = std::make_shared<Field>(
            "MATERIAL_STATE_NEW", FieldDomain::ELEMENT_MP, mp_offset, 1);

        material_state_old->set_zero();
        material_state_new->set_zero();
    }
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
                               bool               fill_nan) {
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

/**
 * Projects element-nodal values onto unique global nodes by weighted averaging.
 *
 * The compiled element-nodal prefix offsets locate the local rows belonging to
 * each dense element. For every connected global node and component, the method
 * accumulates the element-local value multiplied by the scalar weight of its
 * element. The accumulated field is then divided by the sum of participating
 * weights at that node.
 *
 * Elements with zero weight do not participate. Nodes without a contributing
 * element retain the initialized zero row. The operation assumes that every
 * represented element owns exactly one element-nodal row per connectivity node.
 *
 * @param element_nodal Source field in the `ELEMENT_NODAL` domain.
 * @param element_weights Scalar `ELEMENT` field controlling contribution weight.
 * @param name Non-empty name assigned to the returned nodal field.
 * @return Weighted global `NODE` field with the source component count.
 */
Field ModelData::element_nodal_to_nodal(const Field&       element_nodal,
                                        const Field&       element_weights,
                                        const std::string& name) const {
    // Validate field domains, component conventions and compiled prerequisites
    logging::error(element_nodal.domain == FieldDomain::ELEMENT_NODAL,
        "ModelData: element_nodal_to_nodal requires ELEMENT_NODAL input field '",
        element_nodal.name, "'");
    logging::error(element_weights.domain == FieldDomain::ELEMENT,
        "ModelData: element_nodal_to_nodal requires ELEMENT weight field '",
        element_weights.name, "'");
    logging::error(element_weights.components == 1,
        "ModelData: element weight field '", element_weights.name,
        "' must have exactly one component");
    logging::error(element_nodal_offsets != nullptr,
        "ModelData: element nodal offsets are not initialized");
    logging::error(positions != nullptr,
        "ModelData: POSITION field is not initialized");
    logging::error(!name.empty(),
        "ModelData: projected nodal field name cannot be empty");

    // Resolve dense assembly dimensions and validate the complete source layout
    const Index element_count = static_cast<Index>(elements.size());
    const Index node_count    = positions->rows;

    logging::error(element_weights.rows == element_count,
        "ModelData: element weight field '", element_weights.name,
        "' has ", element_weights.rows, " rows, expected ", element_count);

    const Field& offsets      = *element_nodal_offsets;
    const Index expected_rows = static_cast<Index>(offsets(element_count, 0));
    logging::error(element_nodal.rows == expected_rows,
        "ModelData: ELEMENT_NODAL field '", element_nodal.name,
        "' has ", element_nodal.rows, " rows, expected ", expected_rows);

    // Initialize nodal value and weight accumulators for the projection
    Field nodal{name, FieldDomain::NODE, node_count, element_nodal.components};
    nodal.set_zero();
    std::vector<Precision> weight_sum(static_cast<std::size_t>(node_count), Precision(0));

    // Accumulate every participating element-local row at its connected global node
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        const auto& element = elements[static_cast<std::size_t>(elem_idx)];
        if (!element) {
            continue;
        }

        const ID elem_id = element->elem_id;
        logging::error(elem_id >= 0 && static_cast<Index>(elem_id) < element_count,
            "ModelData: element id out of range in element_nodal_to_nodal: ", elem_id);

        const Precision weight = element_weights(static_cast<Index>(elem_id), 0);
        if (weight == Precision(0)) {
            continue;
        }

        const Index offset      = static_cast<Index>(offsets(static_cast<Index>(elem_id), 0));
        const Index next_offset = static_cast<Index>(offsets(static_cast<Index>(elem_id) + 1, 0));
        logging::error(next_offset >= offset,
            "ModelData: invalid element nodal offsets for element ", elem_id);
        logging::error(next_offset - offset == static_cast<Index>(element->n_nodes()),
            "ModelData: element nodal offset span does not match node count for element ", elem_id);

        for (Index local_node = 0; local_node < static_cast<Index>(element->n_nodes()); ++local_node) {
            const Index element_row = offset + local_node;
            const Index node_id     = static_cast<Index>(element->nodes()[local_node]);
            logging::error(node_id < node_count,
                "ModelData: node id out of range in element_nodal_to_nodal: ", node_id);

            for (Index component = 0; component < element_nodal.components; ++component) {
                nodal(node_id, component) += weight * element_nodal(element_row, component);
            }
            weight_sum[static_cast<std::size_t>(node_id)] += weight;
        }
    }

    // Normalize accumulated values by the total element weight at each node
    for (Index node = 0; node < node_count; ++node) {
        const Precision weight = weight_sum[static_cast<std::size_t>(node)];
        if (weight == Precision(0)) {
            continue;
        }
        for (Index component = 0; component < element_nodal.components; ++component) {
            nodal(node, component) /= weight;
        }
    }

    return nodal;
}

} // namespace fem::model
