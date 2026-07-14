#include "writer_res.h"

#include "../model/element/element.h"
#include "../model/model_data.h"

#include <cstddef>
#include <iomanip>

namespace fem {
namespace reader {

namespace {

const char* field_domain_to_string(model::FieldDomain domain) {
    switch (domain) {
        case model::FieldDomain::UNKNOWN:
            return "UNKNOWN";
        case model::FieldDomain::NODE:
            return "NODE";
        case model::FieldDomain::ELEMENT:
            return "ELEMENT";
        case model::FieldDomain::ELEMENT_NODAL:
            return "ELEMENT_NODAL";
        case model::FieldDomain::ELEMENT_IP:
            return "ELEMENT_IP";
    }

    return "UNKNOWN";
}

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

void write_field_values(std::ofstream& file_path,
                        const model::Field& field,
                        Index row) {
    for (Index j = 0; j < field.components; ++j) {
        file_path << std::setw(16) << field(row, j);
    }

    file_path << '\n';
}

void write_dense_field(std::ofstream& file_path,
                       const model::Field& field,
                       const std::string& field_name) {
    const int rows = static_cast<int>(field.rows);
    const int cols = static_cast<int>(field.components);

    write_dense_field_header(file_path,
                             field_name,
                             field.domain,
                             cols,
                             rows);

    file_path.setf(std::ios::scientific, std::ios::floatfield);

    for (Index i = 0; i < field.rows; ++i) {
        write_field_values(file_path, field, i);
    }

    file_path << "END FIELD\n";
}

void write_element_location_field(std::ofstream& file_path,
                                  const model::Field& field,
                                  const std::string& field_name,
                                  const model::ModelData& model_data,
                                  const model::Field& offsets,
                                  const char* offset_name) {
    logging::error(offsets.rows == static_cast<Index>(model_data.max_elems + 1),
                   "ResWriter: ", offset_name, " must have max_elems + 1 rows");

    logging::error(field.rows == static_cast<Index>(offsets(static_cast<Index>(model_data.max_elems), 0)),
                   "ResWriter: field '", field_name, "' row count does not match ", offset_name);

    write_indexed_field_header(file_path,
                               field_name,
                               field.domain,
                               2,
                               static_cast<int>(field.components),
                               static_cast<int>(field.rows));

    file_path.setf(std::ios::scientific, std::ios::floatfield);

    for (Index elem_row = 0; elem_row < static_cast<Index>(model_data.max_elems); ++elem_row) {
        const auto& element = model_data.elements[static_cast<std::size_t>(elem_row)];

        if (!element) {
            continue;
        }

        const Index begin = static_cast<Index>(offsets(elem_row, 0));
        const Index end   = static_cast<Index>(offsets(elem_row + 1, 0));

        logging::error(begin <= end && end <= field.rows,
                       "ResWriter: invalid ", offset_name, " span for element ", elem_row);

        const ID elem_id = element->elem_id;

        for (Index local = 0; local < end - begin; ++local) {
            const Index row = begin + local;

            file_path << std::setw(16) << elem_id
                      << std::setw(16) << local;

            write_field_values(file_path, field, row);
        }
    }

    file_path << "END FIELD\n";
}

} // namespace

ResWriter::ResWriter(const std::string& filename) {
    if (!filename.empty()) {
        open(filename);
    }
}

ResWriter::~ResWriter() {
    close();
}

ResWriter::ResWriter(ResWriter&& other) noexcept
    : file_path(std::move(other.file_path)) {}

ResWriter& ResWriter::operator=(ResWriter&& other) noexcept {
    if (this != &other) {
        close();
        file_path = std::move(other.file_path);
    }

    return *this;
}

void ResWriter::open(const std::string& filename) {
    close();

    file_path.open(filename, std::ios::out | std::ios::trunc);

    logging::error(file_path.is_open(),
                   "ResWriter: failed to open file: ", filename);
}

void ResWriter::close() {
    if (file_path.is_open()) {
        file_path.close();
    }
}

void ResWriter::add_loadcase(int id, WriterStepType step_type) {
    (void) step_type;

    logging::error(file_path.is_open(),
                   "ResWriter: cannot add loadcase: file is not open");

    file_path << "LC " << id << '\n';
}

void ResWriter::write_field(const model::Field& field,
                            const std::string& field_name,
                            const model::ModelData* model_data,
                            Precision frame_value) {
    (void) frame_value;

    logging::error(file_path.is_open(),
                   "ResWriter: cannot write field '", field_name,
                   "': file is not open");

    if (field.domain == model::FieldDomain::ELEMENT_NODAL) {
        logging::error(model_data != nullptr,
                       "ResWriter: ELEMENT_NODAL field '", field_name,
                       "' requires model data");

        logging::error(model_data->element_nodal_offsets != nullptr,
                       "ResWriter: element nodal offsets are not initialized for field '",
                       field_name, "'");

        write_element_location_field(file_path,
                                     field,
                                     field_name,
                                     *model_data,
                                     *model_data->element_nodal_offsets,
                                     "element nodal offsets");
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_IP) {
        logging::error(model_data != nullptr,
                       "ResWriter: ELEMENT_IP field '", field_name,
                       "' requires model data");

        logging::error(model_data->element_ip_offsets != nullptr,
                       "ResWriter: element IP offsets are not initialized for field '",
                       field_name, "'");

        write_element_location_field(file_path,
                                     field,
                                     field_name,
                                     *model_data,
                                     *model_data->element_ip_offsets,
                                     "element IP offsets");
        return;
    }

    write_dense_field(file_path, field, field_name);
}

} // namespace reader
} // namespace fem
