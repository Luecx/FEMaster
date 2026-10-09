/**
 * @file writer_fil.cpp
 * @brief ASCII Abaqus results-file serialization.
 *
 * Records consist of a length word, a key word and typed attributes.
 * Their textual representation is an uninterrupted 80-column stream; logical
 * record boundaries are indicated by '*', not by physical newlines.
 */
#include "writer_fil.h"

#include "../../core/logging.h"
#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cctype>
#include <limits>
#include <string>
#include <system_error>

namespace fem::io::writer {

namespace {

int procedure(WriterStepType type) {
    switch (type) {
        case WriterStepType::Static:         return 1;
        case WriterStepType::Dynamic:        return 11;
        case WriterStepType::Eigenfrequency: return 41;
        case WriterStepType::Buckling:       return 42;
    }
    return 1;
}

std::string normalize(std::string_view name) {
    std::string key;
    for (unsigned char c : name) {
        if (std::isalpha(c)) key.push_back(static_cast<char>(std::toupper(c)));
    }
    return key;
}

} // namespace

FilWriter::FilWriter(const std::string& filename) {
    if (!filename.empty()) open(filename);
}

FilWriter::~FilWriter() {
    close();
}

void FilWriter::open(const std::string& filename) {
    close();
    file_.open(filename, std::ios::binary | std::ios::trunc);
    logging::error(file_.is_open(), "FilWriter: cannot open ", filename);
    line_.clear();
    model_data_ = nullptr;
    model_written_ = false;
    frame_open_ = false;
    rotations_ = false;
    step_ = 1;
    increment_ = 0;
    step_type_ = WriterStepType::Static;
    previous_step_time_ = 0.;
    total_time_offset_ = 0.;
}

void FilWriter::close() {
    if (!file_.is_open()) return;
    end_frame();
    if (!line_.empty()) finish_line();
    file_.close();
}

void FilWriter::append(std::string_view value) {
    while (!value.empty()) {
        const std::size_t count = std::min(value.size(), 80 - line_.size());
        line_.append(value.data(), count);
        value.remove_prefix(count);
        if (line_.size() == 80) {
            file_.write(line_.data(), 80);
            file_.put('\n');
            line_.clear();
        }
    }
}

void FilWriter::finish_line() {
    if (line_.empty()) return;
    line_.append(80 - line_.size(), ' ');
    file_.write(line_.data(), 80);
    file_.put('\n');
    line_.clear();
}

void FilWriter::record(int key, std::size_t attributes) {
    append("*");
    integer(static_cast<long long>(attributes + 2));
    integer(key);
}

void FilWriter::integer(long long value) {
    const std::string digits = std::to_string(value);
    logging::error(digits.size() <= 99, "FilWriter: integer is too long");
    char prefix[3] = {'I', static_cast<char>('0' + digits.size() / 10),
                           static_cast<char>('0' + digits.size() % 10)};
    if (digits.size() < 10) prefix[1] = ' ';
    append(std::string_view(prefix, 3));
    append(digits);
}

void FilWriter::real(double value) {
    logging::error(std::isfinite(value), "FilWriter: non-finite result");
    char buffer[64];
    const auto [last, error] = std::to_chars(buffer, buffer + sizeof(buffer),
                                             value, std::chars_format::scientific, 15);
    logging::error(error == std::errc{}, "FilWriter: cannot format floating point value");
    std::string formatted(buffer, last);
    for (char& c : formatted) if (c == 'e' || c == 'E') c = 'D';
    logging::error(formatted.size() <= 22, "FilWriter: real value exceeds D22.15");
    append("D");
    append(std::string(22 - formatted.size(), ' '));
    append(formatted);
}

void FilWriter::string(std::string_view value) {
    logging::error(value.size() <= 8, "FilWriter: A8 word is too long");
    append("A");
    append(value);
    append(std::string(8 - value.size(), ' '));
}

void FilWriter::write_model_data(const model::ModelData& model_data) {
    logging::error(file_.is_open(), "FilWriter: file not open");
    if (model_written_) return;
    logging::error(model_data.positions != nullptr, "FilWriter: no compiled positions");
    const auto& positions = *model_data.positions;
    logging::error(positions.components >= 3, "FilWriter: positions need 3 components");

    model_data_ = &model_data;

    std::size_t supported_elements = 0;
    for (const auto& element : model_data.elements) {
        if (element && !element->type_name().empty() && element->type_name().size() <= 8)
            ++supported_elements;
    }

    // Model header: release, date (2 words), time, element count, node count,
    // characteristic length. Global dense identifiers are 1-based in .fil.
    record(1921, 7);
    string("FEMASTER");
    string(""); string(""); string("");
    integer(static_cast<long long>(supported_elements));
    integer(static_cast<long long>(positions.rows));
    real(0.0);

    for (Index i = 0; i < positions.rows; ++i) {
        record(1901, 4);
        integer(static_cast<long long>(i) + 1);
        for (int c = 0; c < 3; ++c) real(positions(i, c));
    }

    bool rotations = false;
    for (const auto& element : model_data.elements) {
        if (!element) continue;
        const std::string type = element->type_name();
        // Only types with an Abaqus-compatible A8 name can be represented.
        if (type.empty() || type.size() > 8) continue;
        const Index count = element->n_nodes();
        record(1900, static_cast<std::size_t>(count) + 2);
        integer(static_cast<long long>(element->elem_id) + 1);
        string(type);
        for (Index j = 0; j < count; ++j) {
            integer(static_cast<long long>(element->nodes()[j]) + 1);
        }
        rotations = rotations || type.front() == 'B' || type.front() == 'S';
    }

    // Global DOF position map. Nonrotational models expose three components.
    rotations_ = rotations;
    record(1902, 6);
    for (int i = 0; i < 6; ++i) integer(i < 3 || rotations_ ? i + 1 : 0);
    model_written_ = true;
}

void FilWriter::add_loadcase(int id, WriterStepType type) {
    end_frame();
    if (increment_ != 0) total_time_offset_ += previous_step_time_;
    step_ = id;
    step_type_ = type;
    increment_ = 0;
    previous_step_time_ = 0.;
}

void FilWriter::begin_frame(Precision value) {
    logging::error(file_.is_open(), "FilWriter: file not open");
    logging::error(model_written_, "FilWriter: mesh must be written before results");
    logging::error(step_type_ == WriterStepType::Static,
                   "FilWriter: only static results are supported; harmonic and modal "
                   "steps need additional Abaqus .fil records");
    end_frame();

    const double time = std::isfinite(value) ? static_cast<double>(value) : 0.;
    const double delta = time - previous_step_time_;
    previous_step_time_ = time;
    ++increment_;
    frame_open_ = true;

    record(2000, 21);
    real(total_time_offset_ + time);
    real(time);
    real(0.);
    real(0.);
    integer(procedure(step_type_));
    integer(step_);
    integer(increment_);
    integer(step_type_ == WriterStepType::Eigenfrequency ||
            step_type_ == WriterStepType::Buckling ? 1 : 0);
    real(0.);
    real(0.);
    real(delta);
    for (int i = 0; i < 10; ++i) string("");
}

void FilWriter::end_frame() {
    if (!frame_open_) return;
    record(2001, 0);
    finish_line();  // 2001 must terminate the current 80-character record.
    file_.write(std::string(80, ' ').data(), 80);
    file_.put('\n');
    frame_open_ = false;
}

void FilWriter::output_request(bool nodal, std::string_view type) {
    record(1911, nodal ? 2 : 3);
    integer(nodal ? 1 : 0);
    string("");
    if (!nodal) string(type);
}

void FilWriter::element_header(int element, int point, int location) {
    record(1, 9);
    integer(element);
    integer(point);
    integer(0);        // Solid element: no section point
    integer(location); // 0 = integration point, 2 = element node
    string("");        // No rebar
    integer(3);        // NDI
    integer(3);        // NSHR
    integer(0);
    integer(0);
}

void FilWriter::write_nodal(const model::Field& field, int key) {
    logging::error(field.rows == model_data_->positions->rows,
                   "FilWriter: nodal field row count mismatch");
    // Abaqus 1902 maps active DOFs: do not emit artificial zero rotations
    // from FEMaster's six-component storage for a solid-only model.
    const Index count = key == 201 ? 1 : std::min<Index>(field.components, rotations_ ? 6 : 3);
    logging::error(field.components >= count, "FilWriter: insufficient nodal components");
    output_request(true);
    for (Index row = 0; row < field.rows; ++row) {
        record(key, static_cast<std::size_t>(count) + 1);
        integer(static_cast<long long>(row) + 1);
        for (Index c = 0; c < count; ++c) real(field(row, c));
    }
}

void FilWriter::write_element(const model::Field& field, int key) {
    const bool nodal = field.domain == model::FieldDomain::ELEMENT_NODAL;
    for (const auto& element : model_data_->elements) {
        if (!element) continue;
        const std::string type = element->type_name();
        // A six-component global Cartesian tensor is not a valid generic
        // section-point representation for shells, beams or trusses.
        if (type.compare(0, 3, "C3D") != 0) continue;
        const Index count = nodal ? element->n_nodes() : element->num_ip();
        const Index offset = nodal ? element->elem_nodal_offset : element->elem_ip_offset;
        output_request(false, type);
        for (Index i = 0; i < count; ++i) {
            const Index row = offset + i;
            logging::error(row < field.rows, "FilWriter: element field row out of range");
            element_header(static_cast<int>(element->elem_id) + 1,
                           nodal ? static_cast<int>(element->nodes()[i]) + 1
                                 : static_cast<int>(i) + 1,
                           nodal ? 2 : 0);
            record(key, 6);
            // FEMaster: XX,YY,ZZ,YZ,ZX,XY; Abaqus: 11,22,33,12,13,23.
            for (int c : {0, 1, 2, 5, 4, 3}) real(field(row, c));
        }
    }
}

void FilWriter::write_field(const model::Field& field, const std::string& name,
                            const model::ModelData* data, Precision value) {
    if (!file_.is_open()) return;
    if (!model_written_) {
        logging::error(data != nullptr, "FilWriter: field requires compiled model data");
        write_model_data(*data);
    }

    const std::string key = normalize(name);
    int record_key = 0;
    if      (key == "DISPLACEMENT")       record_key = 101;
    else if (key == "VELOCITY")           record_key = 102;
    else if (key == "ACCELERATION")       record_key = 103;
    else if (key == "REACTIONFORCES")     record_key = 104;
    else if (key == "EXTERNALFORCES")     record_key = 106;
    else if (key == "TEMPERATURE")        record_key = 201;
    else if (key == "STRESS")             record_key = 11;
    else if (key == "STRAIN")             record_key = 21;

    const bool nodal = field.domain == model::FieldDomain::NODE && record_key >= 101;
    const bool tensor = (field.domain == model::FieldDomain::ELEMENT_NODAL ||
                         field.domain == model::FieldDomain::ELEMENT_IP) &&
                         (record_key == 11 || record_key == 21) &&
                         field.components == 6;
    if (!nodal && !tensor) return;

    if (!frame_open_) begin_frame(value);
    if (nodal) write_nodal(field, record_key);
    else       write_element(field, record_key);
}

} // namespace fem::io::writer
