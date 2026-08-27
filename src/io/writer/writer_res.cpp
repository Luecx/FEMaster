/**
 * @file writer_res.cpp
 * @brief Implements the native FEMaster result writer.
 *
 * NODE and ELEMENT fields are written as explicitly indexed fields. The writer
 * converts dense assembly identifiers to semantic result identifiers once from
 * the compiled Part/Instance topology and caches the resulting strings.
 *
 * Numerical field data is formatted into thread-local text buffers before being
 * written to disk. Large fields may therefore be formatted in parallel while
 * preserving their deterministic row order in the final RES file.
 *
 * The implicit default Instance preserves the original part-local identifier:
 *
 *     17
 *
 * Explicit Instances use a qualified identifier:
 *
 *     bolt.17
 *
 * The same element identifiers are reused for ELEMENT, ELEMENT_NODAL,
 * ELEMENT_IP and ELEMENT_MP fields.
 *
 * @see ResWriter
 * @see model::ModelData
 *
 * @author Finn Eggers
 * @date 20.08.2026
 */

#include "writer_res.h"

#include "../../core/config.h"
#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace fem {
namespace io {
namespace writer {

namespace {

/**
 * Appends a right-aligned string value to a RES output buffer.
 *
 * The requested width is a minimum width and therefore follows the behavior of
 * `std::setw`: values longer than the requested width are written completely.
 */
void append_string(std::string& buffer,
                   std::string_view value,
                   std::size_t width = 16) {
    if (value.size() < width) {
        buffer.append(width - value.size(), ' ');
    }

    buffer.append(value);
}

/**
 * Appends a right-aligned integer value to a RES output buffer.
 *
 * Integer conversion uses `std::to_chars` to avoid locale handling and stream
 * formatting overhead in the inner result-writing loops.
 */
template<typename T>
void append_integer(std::string& buffer,
                    T value,
                    std::size_t width = 16) {
    char tmp[32];

    const auto [end, ec] = std::to_chars(tmp, tmp + sizeof(tmp), value);

    logging::error(ec == std::errc{},
        "ResWriter: failed to format integer");

    const std::size_t length = static_cast<std::size_t>(end - tmp);

    if (length < width) {
        buffer.append(width - length, ' ');
    }

    buffer.append(tmp, length);
}

/**
 * Appends a right-aligned floating-point value in scientific notation.
 *
 * Six digits after the decimal point preserve the formatting previously
 * produced by an `std::ofstream` in scientific mode with its default precision.
 */
void append_precision(std::string& buffer,
                      Precision value,
                      std::size_t width = 16) {
    char tmp[64];

    const auto [end, ec] = std::to_chars(
        tmp,
        tmp + sizeof(tmp),
        value,
        std::chars_format::scientific,
        6
    );

    logging::error(ec == std::errc{},
        "ResWriter: failed to format floating-point value");

    const std::size_t length = static_cast<std::size_t>(end - tmp);

    if (length < width) {
        buffer.append(width - length, ' ');
    }

    buffer.append(tmp, length);
}

/**
 * Appends all numerical components of one field row and terminates the row.
 */
void append_field_values(std::string& buffer,
                         const model::Field& field,
                         Index row) {
    for (Index component = 0; component < field.components; ++component) {
        append_precision(buffer, field(row, component));
    }

    buffer.push_back('\n');
}

/**
 * Selects the number of workers used to format one RES field.
 *
 * Small fields remain single-threaded to avoid OpenMP and buffer-management
 * overhead. Large fields use at most one worker per 4096 work items and never
 * exceed either the configured FEMaster thread count or 32 writer threads.
 */
ID output_thread_count(Index work_items) {
#ifndef _OPENMP
    (void) work_items;
    return 1;
#else
    const ID useful_threads = static_cast<ID>(
        std::min<Index>(std::max<Index>(1, work_items / 4096), 32)
    );

    return std::max<ID>(1, std::min<ID>(global_config.max_threads, useful_threads));
#endif
}

/**
 * Creates one independent output buffer for every formatting worker.
 */
std::vector<std::string> create_thread_buffers(ID num_threads) {
    std::vector<std::string> buffers(static_cast<std::size_t>(num_threads));

    for (auto& buffer : buffers) {
        buffer.reserve(16 * 1024 * 1024);
    }

    return buffers;
}

/**
 * Writes completed worker buffers sequentially to preserve field row order.
 */
void write_thread_buffers(std::ofstream& file_path,
                          const std::vector<std::string>& buffers) {
    for (const auto& buffer : buffers) {
        file_path.write(
            buffer.data(),
            static_cast<std::streamsize>(buffer.size())
        );
    }
}

/**
 * Converts a field domain to its textual RES representation.
 */
const char* field_domain_to_string(model::FieldDomain domain) {
    switch (domain) {
        case model::FieldDomain::UNKNOWN:       return "UNKNOWN";
        case model::FieldDomain::NODE:          return "NODE";
        case model::FieldDomain::ELEMENT:       return "ELEMENT";
        case model::FieldDomain::ELEMENT_NODAL: return "ELEMENT_NODAL";
        case model::FieldDomain::ELEMENT_IP:    return "ELEMENT_IP";
        case model::FieldDomain::ELEMENT_MP:    return "ELEMENT_MP";
    }

    return "UNKNOWN";
}

/**
 * Writes the header of a dense RES field.
 */
void write_dense_field_header(std::ofstream& file_path,
                              const std::string& field_name,
                              model::FieldDomain domain,
                              int cols,
                              int rows) {
    file_path << "FIELD, NAME=" << field_name
              << ", TYPE=" << field_domain_to_string(domain)
              << ", COLS=" << cols
              << ", ROWS=" << rows << "\n";
}

/**
 * Writes the header of an explicitly indexed RES field.
 */
void write_indexed_field_header(std::ofstream& file_path,
                                const std::string& field_name,
                                model::FieldDomain domain,
                                int index_cols,
                                int value_cols,
                                int rows) {
    file_path << "FIELD, NAME=" << field_name
              << ", TYPE=" << field_domain_to_string(domain)
              << ", INDEX_COLS=" << index_cols
              << ", VALUE_COLS=" << value_cols
              << ", ROWS=" << rows << "\n";
}

/**
 * Writes a field without explicit entity identifiers.
 *
 * Numerical rows are divided into contiguous ranges. Each worker formats one
 * range into an independent string buffer, after which the buffers are written
 * sequentially so the original dense row ordering is retained.
 *
 * This path is retained for UNKNOWN fields. All entity-based NODE and ELEMENT
 * domains are written through indexed output instead.
 */
void write_dense_field(std::ofstream& file_path,
                       const model::Field& field,
                       const std::string& field_name) {
    write_dense_field_header(
        file_path,
        field_name,
        field.domain,
        static_cast<int>(field.components),
        static_cast<int>(field.rows)
    );

    // Allocate one contiguous output buffer for every formatting worker
    const ID num_threads = output_thread_count(field.rows);
    auto thread_buffers = create_thread_buffers(num_threads);

    // Format contiguous row ranges independently
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
#endif
    for (ID thread = 0; thread < num_threads; ++thread) {
        auto& buffer = thread_buffers[static_cast<std::size_t>(thread)];

        const Index begin = field.rows * thread       / num_threads;
        const Index end   = field.rows * (thread + 1) / num_threads;

        for (Index row = begin; row < end; ++row) {
            append_field_values(buffer, field, row);
        }
    }

    // Preserve dense row order while committing the formatted text to disk
    write_thread_buffers(file_path, thread_buffers);

    file_path << "END FIELD\n";
}

/**
 * Writes one NODE or ELEMENT field with an explicit semantic entity identifier.
 *
 * Field rows remain in dense assembly order. Each worker receives one contiguous
 * row range and formats both the semantic identifier and numerical field values
 * into its private output buffer. Worker buffers are subsequently written in
 * their original range order.
 *
 * @param file_path Active RES output stream.
 * @param field Field whose rows are written.
 * @param field_name Name stored in the RES field header.
 * @param ids Semantic identifier corresponding to every dense field row.
 */
void write_indexed_field(std::ofstream& file_path,
                         const model::Field& field,
                         const std::string& field_name,
                         const std::vector<std::string>& ids) {
    logging::error(ids.size() == static_cast<std::size_t>(field.rows),
        "ResWriter: id mapping does not match rows of field '", field_name, "'");

    write_indexed_field_header(
        file_path,
        field_name,
        field.domain,
        1,
        static_cast<int>(field.components),
        static_cast<int>(field.rows)
    );

    // Allocate one contiguous output buffer for every formatting worker
    const ID num_threads = output_thread_count(field.rows);
    auto thread_buffers = create_thread_buffers(num_threads);

    // Format semantic identifiers and field rows independently
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
#endif
    for (ID thread = 0; thread < num_threads; ++thread) {
        auto& buffer = thread_buffers[static_cast<std::size_t>(thread)];

        const Index begin = field.rows * thread       / num_threads;
        const Index end   = field.rows * (thread + 1) / num_threads;

        for (Index row = begin; row < end; ++row) {
            append_string(buffer, ids[static_cast<std::size_t>(row)]);
            append_field_values(buffer, field, row);
        }
    }

    // Preserve dense row order while committing the formatted text to disk
    write_thread_buffers(file_path, thread_buffers);

    file_path << "END FIELD\n";
}

/**
 * Writes an element-location field with semantic element identifiers.
 *
 * The first index identifies the element using its RES identifier. The second
 * index retains the existing local node, integration-point or material-point
 * enumeration within that element.
 *
 * Elements are split into contiguous ranges so every worker owns complete
 * elements and their associated local rows. This preserves the established
 * element-major ordering while allowing numerical text conversion to run in
 * parallel without synchronization.
 *
 * @param file_path Active RES output stream.
 * @param field Element-location field whose rows are written.
 * @param field_name Name stored in the RES field header.
 * @param model_data Compiled model providing the element enumeration.
 * @param element_ids Semantic identifier corresponding to every dense element.
 * @param offsets Prefix offsets mapping elements to field-row ranges.
 * @param offset_name Human-readable offset name used in diagnostics.
 */
void write_element_location_field(std::ofstream& file_path,
                                  const model::Field& field,
                                  const std::string& field_name,
                                  const model::ModelData& model_data,
                                  const std::vector<std::string>& element_ids,
                                  const model::Field& offsets,
                                  const char* offset_name) {
    const Index element_count = static_cast<Index>(model_data.elements.size());

    logging::error(offsets.rows == element_count + 1,
        "ResWriter: ", offset_name, " must have element_count + 1 rows");
    logging::error(field.rows == static_cast<Index>(offsets(element_count, 0)),
        "ResWriter: field '", field_name, "' row count does not match ", offset_name);
    logging::error(element_ids.size() == model_data.elements.size(),
        "ResWriter: element id mapping does not match compiled elements");

    write_indexed_field_header(
        file_path,
        field_name,
        field.domain,
        2,
        static_cast<int>(field.components),
        static_cast<int>(field.rows)
    );

    // Keep complete element ranges together inside one formatting worker
    const ID num_threads = output_thread_count(element_count);
    auto thread_buffers = create_thread_buffers(num_threads);

    // Format each contiguous element range into its private output buffer
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
#endif
    for (ID thread = 0; thread < num_threads; ++thread) {
        auto& buffer = thread_buffers[static_cast<std::size_t>(thread)];

        const Index start_index = element_count * thread       / num_threads;
        const Index end_index   = element_count * (thread + 1) / num_threads;

        for (Index elem_row = start_index; elem_row < end_index; ++elem_row) {
            const auto& element =
                model_data.elements[static_cast<std::size_t>(elem_row)];

            if (!element) {
                continue;
            }

            const Index begin = static_cast<Index>(offsets(elem_row, 0));
            const Index end   = static_cast<Index>(offsets(elem_row + 1, 0));

            logging::error(begin <= end && end <= field.rows,
                "ResWriter: invalid ", offset_name, " span for element ", elem_row);

            const std::string& elem_id =
                element_ids[static_cast<std::size_t>(elem_row)];

            for (Index local = 0; local < end - begin; ++local) {
                const Index row = begin + local;

                append_string(buffer, elem_id);
                append_integer(buffer, local);
                append_field_values(buffer, field, row);
            }
        }
    }

    // Preserve element-major ordering while committing the buffers to disk
    write_thread_buffers(file_path, thread_buffers);

    file_path << "END FIELD\n";
}

} // namespace

/**
 * Constructs a RES writer and optionally opens the target file.
 */
ResWriter::ResWriter(const std::string& filename) {
    if (!filename.empty()) {
        open(filename);
    }
}

/**
 * Closes the active result file.
 */
ResWriter::~ResWriter() {
    close();
}

/**
 * Moves an active writer together with its cached semantic entity identifiers.
 */
ResWriter::ResWriter(ResWriter&& other) noexcept
    : file_path      (std::move(other.file_path)),
      node_ids       (std::move(other.node_ids)),
      element_ids    (std::move(other.element_ids)),
      ids_initialized(other.ids_initialized) {
    other.ids_initialized = false;
}

/**
 * Moves an active writer together with its cached semantic entity identifiers.
 */
ResWriter& ResWriter::operator=(ResWriter&& other) noexcept {
    if (this != &other) {
        close();

        file_path       = std::move(other.file_path);
        node_ids        = std::move(other.node_ids);
        element_ids     = std::move(other.element_ids);
        ids_initialized = other.ids_initialized;

        other.ids_initialized = false;
    }

    return *this;
}

/**
 * Opens a new RES file and clears all writer-local identifier mappings.
 */
void ResWriter::open(const std::string& filename) {
    close();

    file_path.open(filename, std::ios::out | std::ios::trunc);

    logging::error(file_path.is_open(),
        "ResWriter: failed to open file: ", filename);

    node_ids.clear();
    element_ids.clear();
    ids_initialized = false;
}

/**
 * Closes the active RES file.
 */
void ResWriter::close() {
    if (file_path.is_open()) {
        file_path.close();
    }
}

/**
 * Starts a new result loadcase.
 */
void ResWriter::add_loadcase(int id, WriterStepType step_type) {
    (void) step_type;

    logging::error(file_path.is_open(),
        "ResWriter: cannot add loadcase: file is not open");

    file_path << "LC " << id << '\n';
}

/**
 * Builds the semantic RES identifiers for compiled nodes and elements.
 *
 * Nodes are resolved directly through `ModelData::node_mapping`, which maps a
 * dense assembly node to its owning Instance and original part-local node id.
 *
 * Elements are resolved by reversing each compiled Instance's existing
 * `element_ids` map from local to dense identifiers.
 *
 * The default Instance retains bare local identifiers for compatibility with
 * existing FEMaster result files. Explicit Instances are qualified as
 * `instance_name.local_id`.
 */
void ResWriter::initialize_ids(const model::ModelData& model_data) {
    logging::error(!ids_initialized,
        "ResWriter: result ids have already been initialized");

    // Build dense node -> semantic node identifier
    node_ids.resize(model_data.node_mapping.size());

    for (std::size_t row = 0; row < model_data.node_mapping.size(); ++row) {
        const auto& [instance, local_id] = model_data.node_mapping[row];

        logging::error(instance != nullptr,
            "ResWriter: dense node ", row, " has no source instance");

        node_ids[row] = instance->instance_id == 0
            ? std::to_string(local_id)
            : instance->name + "." + std::to_string(local_id);
    }

    // Reverse every Instance's local -> dense element map
    element_ids.assign(model_data.elements.size(), "");

    for (const auto& [instance_name, instance] : model_data.instances) {
        logging::error(instance != nullptr,
            "ResWriter: instance ", instance_name, " is null");

        for (const auto& [local_id, global_id] : instance->element_ids) {
            const std::size_t row = static_cast<std::size_t>(global_id);

            logging::error(row < element_ids.size(),
                "ResWriter: element ", local_id, " in instance ", instance_name,
                " maps outside compiled element storage");
            logging::error(element_ids[row].empty(),
                "ResWriter: dense element ", global_id,
                " has multiple semantic identifiers");

            element_ids[row] = instance->instance_id == 0
                ? std::to_string(local_id)
                : instance_name + "." + std::to_string(local_id);
        }
    }

    // Every compiled element must have a semantic source identifier
    for (std::size_t row = 0; row < model_data.elements.size(); ++row) {
        if (!model_data.elements[row]) {
            continue;
        }

        logging::error(!element_ids[row].empty(),
            "ResWriter: dense element ", row,
            " has no semantic element identifier");
    }

    ids_initialized = true;
}

/**
 * Writes one FEMaster result field.
 *
 * All entity-based field domains are explicitly indexed. The semantic node and
 * element mappings are initialized lazily on the first such write and retained
 * for subsequent fields.
 *
 * ELEMENT_NODAL, ELEMENT_IP and ELEMENT_MP fields still require ModelData on
 * every write because their dense row spans are defined by the corresponding
 * element-offset fields.
 */
void ResWriter::write_field(const model::Field& field,
                            const std::string& field_name,
                            const model::ModelData* model_data,
                            Precision frame_value) {
    (void) frame_value;

    logging::error(file_path.is_open(),
        "ResWriter: cannot write field '", field_name,
        "': file is not open");

    const bool indexed =
           field.domain == model::FieldDomain::NODE
        || field.domain == model::FieldDomain::ELEMENT
        || field.domain == model::FieldDomain::ELEMENT_NODAL
        || field.domain == model::FieldDomain::ELEMENT_IP
        || field.domain == model::FieldDomain::ELEMENT_MP;

    if (indexed && !ids_initialized) {
        logging::error(model_data != nullptr,
            "ResWriter: field '", field_name,
            "' requires model data to initialize result ids");

        initialize_ids(*model_data);
    }

    if (field.domain == model::FieldDomain::NODE) {
        write_indexed_field(file_path, field, field_name, node_ids);
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT) {
        write_indexed_field(file_path, field, field_name, element_ids);
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_NODAL) {
        logging::error(model_data != nullptr,
            "ResWriter: ELEMENT_NODAL field '", field_name,
            "' requires model data");
        logging::error(model_data->element_nodal_offsets != nullptr,
            "ResWriter: element nodal offsets are not initialized for field '",
            field_name, "'");

        write_element_location_field(
            file_path,
            field,
            field_name,
            *model_data,
            element_ids,
            *model_data->element_nodal_offsets,
            "element nodal offsets"
        );
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_IP) {
        logging::error(model_data != nullptr,
            "ResWriter: ELEMENT_IP field '", field_name,
            "' requires model data");
        logging::error(model_data->element_ip_offsets != nullptr,
            "ResWriter: element IP offsets are not initialized for field '",
            field_name, "'");

        write_element_location_field(
            file_path,
            field,
            field_name,
            *model_data,
            element_ids,
            *model_data->element_ip_offsets,
            "element IP offsets"
        );
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_MP) {
        logging::error(model_data != nullptr,
            "ResWriter: ELEMENT_MP field '", field_name,
            "' requires model data");
        logging::error(model_data->element_mp_offsets != nullptr,
            "ResWriter: element MP offsets are not initialized for field '",
            field_name, "'");

        write_element_location_field(
            file_path,
            field,
            field_name,
            *model_data,
            element_ids,
            *model_data->element_mp_offsets,
            "element MP offsets"
        );
        return;
    }

    write_dense_field(file_path, field, field_name);
}

} // namespace writer
} // namespace io
} // namespace fem