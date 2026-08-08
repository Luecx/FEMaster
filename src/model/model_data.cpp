/**
 * @file model_data.cpp
 * @brief Implements model field allocation and element-local enumeration.
 *
 * `ModelData` derives dense model-field row counts from the active model
 * domains and enumerates element-nodal, integration-point and material-point
 * offsets. Material-point enumeration also creates the always-bound default
 * `MATERIAL_STATE` field used by direct element-to-material state addressing;
 * nonlinear analyses may temporarily replace this binding with their active
 * trial buffer.
 *
 * @see ModelData
 * @see Field
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "model_data.h"
#include "element/element.h"

namespace fem {
namespace model {

/**
 * Returns the dense row count associated with one field domain.
 *
 * Node and element domains use their fixed model capacities. Flattened element-
 * nodal, integration-point and material-point domains use the sentinel entry of
 * their prefix-offset field and are therefore unavailable before topology
 * enumeration. Unknown domains and empty offset-dependent domains are rejected.
 *
 * @param domain Semantic field domain whose row count is requested.
 * @return Number of rows required by a field in that domain.
 */
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
        case FieldDomain::ELEMENT_MP:
            logging::error(element_mp_offsets != nullptr,
                           "ModelData: element MP offset field is not initialized");
            logging::error((*element_mp_offsets)(static_cast<Index>(max_elems)) > 0,
                           "ModelData: no material points are available, cannot allocate ELEMENT_MP fields");
            return static_cast<Index>((*element_mp_offsets)(static_cast<Index>(max_elems)));
    }
    logging::error(false, "ModelData: unknown field domain");
    return 0;
}

/**
 * Enumerates all element-local nodal, integration-point and material-point rows.
 *
 * Three prefix-offset fields are built in global element-array order. Each
 * existing element receives the first row of its contiguous span, and the final
 * sentinel entry stores the total row count. Material points are ordered first
 * by element, then by integration point and finally by the section-defined
 * material point within that integration point.
 *
 * After enumeration, a one-component default material-state field is created
 * for every material-point row. This guarantees valid direct addressing even
 * for stateless materials. A nonlinear analysis may replace that binding with a
 * wider trial field after querying the assigned constitutive state sizes.
 *
 * The topology must be complete and enumeration may run only once because
 * element offsets become persistent addressing invariants.
 */
void ModelData::initialize_element_enumeration() {
    // Reject repeated enumeration because existing fields may already depend on
    // the established element-local row layout
    logging::error(element_nodal_offsets == nullptr &&
                   element_ip_offsets    == nullptr &&
                   element_mp_offsets    == nullptr,
                   "ModelData: element enumeration has already been initialized");

    // Allocate prefix arrays with one terminal sentinel beyond the element range
    element_nodal_offsets = std::make_shared<Field>(
        "ELEMENT_NODAL_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);
    element_ip_offsets = std::make_shared<Field>(
        "ELEMENT_IP_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);
    element_mp_offsets = std::make_shared<Field>(
        "ELEMENT_MP_OFFSETS", FieldDomain::ELEMENT, static_cast<Index>(max_elems + 1), 1);

    Index nodal_offset = 0;
    Index ip_offset    = 0;
    Index mp_offset    = 0;

    // Assign every element the first row of each flattened local-data span
    for (Index row = 0; row < static_cast<Index>(max_elems); ++row) {
        (*element_nodal_offsets)(row) = static_cast<Precision>(nodal_offset);
        (*element_ip_offsets)(row)    = static_cast<Precision>(ip_offset);
        (*element_mp_offsets)(row)    = static_cast<Precision>(mp_offset);

        const ElementPtr& element = elements[row];
        if (element != nullptr) {
            element->elem_nodal_offset = static_cast<ID>(nodal_offset);
            element->elem_ip_offset    = static_cast<ID>(ip_offset);
            element->elem_mp_offset    = static_cast<ID>(mp_offset);

            nodal_offset += static_cast<Index>(element->n_nodes());
            ip_offset    += static_cast<Index>(element->num_ip());
            mp_offset    += static_cast<Index>(element->num_ip()) * element->num_mp_per_ip();
        }
    }

    // Store total row counts in the sentinel entries and cache the resulting
    // integration-point and material-point capacities
    const Index sentinel = static_cast<Index>(max_elems);
    (*element_nodal_offsets)(sentinel) = static_cast<Precision>(nodal_offset);
    (*element_ip_offsets)(sentinel)    = static_cast<Precision>(ip_offset);
    (*element_mp_offsets)(sentinel)    = static_cast<Precision>(mp_offset);
    max_integration_points = static_cast<ID>(ip_offset);
    max_material_points    = static_cast<ID>(mp_offset);

    // Keep one directly addressable state row for every enumerated material
    // point. Stateless constitutive laws simply ignore the single dummy
    // component. Nonlinear state management replaces this field with the
    // correctly sized active trial buffer when history variables are present.
    if (max_material_points > 0) {
        material_state = std::make_shared<Field>(
            "MATERIAL_STATE",
            FieldDomain::ELEMENT_MP,
            static_cast<Index>(max_material_points),
            1
        );
        material_state->set_zero();
    }
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

/**
 * Creates field storage for one semantic domain and optionally registers it by
 * name.
 *
 * Registered creation is idempotent: an existing field is returned after its
 * domain and component count have been validated. Unregistered creation always
 * produces independent temporary storage. New fields derive their row count
 * from `field_rows()` and may be initialized with NaN to expose missing writes.
 *
 * @param name Non-empty field name.
 * @param domain Semantic domain controlling the row count.
 * @param components Number of scalar components stored in every row.
 * @param fill_nan Initialize new storage with NaN when `true`.
 * @param reg Register and reuse the field in the model field dictionary.
 * @return Shared field storage.
 */
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

/**
 * Projects a flattened element-nodal field onto shared global nodes by weighted
 * averaging.
 *
 * Every element contributes each of its local nodal rows with the scalar weight
 * stored for that element. Contributions are accumulated component-wise in
 * global node order and normalized by the sum of incident non-zero weights.
 * Nodes without a contribution remain zero. Prefix-offset consistency and node
 * identifiers are validated before accessing the flattened input rows.
 *
 * @param element_nodal Source field in `ELEMENT_NODAL` ordering.
 * @param element_weights Scalar weight for every global element slot.
 * @param name Name of the returned nodal field.
 * @return Weighted nodal projection with the source component count.
 */
Field ModelData::element_nodal_to_nodal(const Field& element_nodal,
                                        const Field& element_weights,
                                        const std::string& name) const {
    // Validate domains, dimensions and the element-nodal addressing metadata
    logging::error(element_nodal.domain == FieldDomain::ELEMENT_NODAL,
                   "ModelData: element_nodal_to_nodal requires ELEMENT_NODAL input field '",
                   element_nodal.name, "'");
    logging::error(element_weights.domain == FieldDomain::ELEMENT,
                   "ModelData: element_nodal_to_nodal requires ELEMENT weight field '",
                   element_weights.name, "'");
    logging::error(element_weights.components == 1,
                   "ModelData: element weight field '", element_weights.name,
                   "' must have exactly one component");
    logging::error(element_weights.rows == static_cast<Index>(max_elems),
                   "ModelData: element weight field '", element_weights.name,
                   "' has ", element_weights.rows, " rows, expected ", max_elems);
    logging::error(element_nodal_offsets != nullptr,
                   "ModelData: element nodal offset field is not initialized");
    logging::error(!name.empty(), "ModelData: projected nodal field name cannot be empty");

    const Field& offsets = *element_nodal_offsets;
    const Index expected_rows = static_cast<Index>(offsets(static_cast<Index>(max_elems), 0));
    logging::error(element_nodal.rows == expected_rows,
                   "ModelData: ELEMENT_NODAL field '", element_nodal.name,
                   "' has ", element_nodal.rows, " rows, expected ", expected_rows);

    // Initialize global accumulation and one scalar normalization weight per node
    Field nodal{name, FieldDomain::NODE, static_cast<Index>(max_nodes), element_nodal.components};
    nodal.set_zero();
    std::vector<Precision> weight_sum(static_cast<std::size_t>(max_nodes), Precision(0));

    // Accumulate every active element-local node into its shared global node
    for (Index elem_idx = 0; elem_idx < static_cast<Index>(max_elems); ++elem_idx) {
        const ElementPtr& element = elements[elem_idx];
        if (!element) {
            continue;
        }

        const ID elem_id = element->elem_id;
        logging::error(elem_id >= 0 && elem_id < max_elems,
                       "ModelData: element id out of range in element_nodal_to_nodal: ", elem_id);

        const Precision weight = element_weights(static_cast<Index>(elem_id), 0);
        if (weight == Precision(0)) {
            continue;
        }

        const Index offset = static_cast<Index>(offsets(static_cast<Index>(elem_id), 0));
        const Index next_offset = static_cast<Index>(offsets(static_cast<Index>(elem_id) + 1, 0));
        logging::error(next_offset >= offset,
                       "ModelData: invalid element nodal offsets for element ", elem_id);
        logging::error(next_offset - offset == static_cast<Index>(element->n_nodes()),
                       "ModelData: element nodal offset span does not match element node count for element ",
                       elem_id);

        for (Index local_node = 0; local_node < static_cast<Index>(element->n_nodes()); ++local_node) {
            const Index element_row = offset + local_node;
            const Index node_id = static_cast<Index>(element->nodes()[local_node]);
            logging::error(node_id < static_cast<Index>(max_nodes),
                           "ModelData: node id out of range in element_nodal_to_nodal: ", node_id);

            for (Index component = 0; component < element_nodal.components; ++component) {
                nodal(node_id, component) += weight * element_nodal(element_row, component);
            }
            weight_sum[static_cast<std::size_t>(node_id)] += weight;
        }
    }

    // Normalize only nodes that received at least one non-zero contribution
    for (Index node = 0; node < static_cast<Index>(max_nodes); ++node) {
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
} // namespace model
} // namespace fem
