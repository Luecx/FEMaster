/**
 * @file writer_fil.cpp
 * @brief Implements the binary Abaqus results-file (.fil) writer.
 *
 * @see writer_fil.h
 *
 * @author Finn Eggers
 * @date 17.11.2026
 */

#include "writer_fil.h"

#include "../../math/pi.h"
#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <ctime>

namespace fem {
namespace io {
namespace writer {

/**
 * Constructs a FIL writer and optionally opens its output file immediately.
 */
FilWriter::FilWriter(const std::string& filename) {
    if (!filename.empty()) {
        open(filename);
    }
}

/**
 * Finalizes the open increment, flushes the last block and closes the file.
 */
FilWriter::~FilWriter() {
    close();
}

/**
 * Opens a new binary FIL file and resets all writer-local state.
 */
void FilWriter::open(const std::string& filename) {
    close();

    file_.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);

    logging::error(file_.is_open(),
        "FilWriter: failed to open file: ", filename);

    node_ids_.clear();
    element_ids_.clear();
    element_types_.clear();
    record_.clear();
    active_dofs_.fill(0);

    block_fill_         = 0;
    model_data_written_ = false;
    increment_open_     = false;
    current_step_       = 1;
    current_frame_      = -1;
    current_step_type_  = WriterStepType::Static;
}

/**
 * Terminates the open increment with a padded 2001 record and closes the file.
 * A file without increments is completed by zero-padding its last block.
 */
void FilWriter::close() {
    if (!file_.is_open()) {
        return;
    }

    if (increment_open_) {
        end_increment();
    }

    flush_block();
    file_.close();
}

/**
 * Selects the step metadata used by subsequent increments. The open increment
 * of the previous step is terminated first.
 */
void FilWriter::add_loadcase(int id, WriterStepType step_type) {
    logging::error(file_.is_open(),
        "FilWriter: cannot add loadcase: file is not open");

    if (increment_open_) {
        end_increment();
    }

    current_step_      = id > 0 ? id : 1;
    current_step_type_ = step_type;
    current_frame_     = -1;
}

/**
 * Builds the persistent external numbering and Abaqus type names and writes
 * the model definition.
 *
 * @param structural_dofs Complete node-by-six structural DOF activation mask.
 */
void FilWriter::write_model_data(const model::ModelData& model_data,
                                 const SystemDofs&       structural_dofs) {
    logging::error(file_.is_open(),
        "FilWriter: cannot write model data: file is not open");
    logging::error(model_data.positions != nullptr,
        "FilWriter: model_data.positions is not initialized");
    logging::error(model_data.node_mapping.size()    == static_cast<std::size_t>(model_data.positions->rows),
        "FilWriter: node mapping does not match positions field");
    logging::error(model_data.element_mapping.size() == model_data.elements.size(),
        "FilWriter: element mapping does not match compiled elements");

    if (model_data_written_) {
        return;
    }

    // Same per-instance label ranges as the FRD writer
    node_ids_.resize(model_data.node_mapping.size());
    for (std::size_t i = 0; i < node_ids_.size(); ++i) {
        const auto& [instance, local_id] = model_data.node_mapping[i];
        logging::error(instance != nullptr, "FilWriter: dense node ", i, " has no source instance");
        node_ids_[i] = std::int64_t(100000000) * instance->instance_id + static_cast<std::int64_t>(local_id);
    }

    element_ids_  .assign(model_data.element_mapping.size(), 0);
    element_types_.assign(model_data.element_mapping.size(), std::string{});

    std::string omitted;
    for (std::size_t i = 0; i < element_ids_.size(); ++i) {
        if (!model_data.elements[i]) continue;
        const auto& [instance, local_id] = model_data.element_mapping[i];
        logging::error(instance != nullptr, "FilWriter: dense element ", i, " has no source instance");
        element_ids_[i] = std::int64_t(100000000) * instance->instance_id + static_cast<std::int64_t>(local_id);

        const std::string name = model_data.elements[i]->type_name();
        element_types_[i]      = abaqus_type(name);
        if (element_types_[i].empty() && omitted.find(" " + name + " ") == std::string::npos) {
            omitted += " " + name + " ";
        }
    }

    logging::warning(omitted.empty(),
        "FilWriter: element types without Abaqus equivalent are omitted from the .fil mesh:", omitted);

    write_header  (model_data);
    write_elements(model_data);
    write_nodes   (model_data);
    write_dofs    (structural_dofs);

    model_data_written_ = true;
}

/**
 * Writes one supported result field into the increment of its frame.
 *
 * Frames are identified by the numeric `_<n>` suffix appended by the output
 * request handler; fields without suffix belong to frame one. A new frame
 * terminates the previous increment and opens the next one.
 */
void FilWriter::write_field(const model::Field&     field,
                            const std::string&      field_name,
                            const model::ModelData* model_data,
                            Precision               frame_value) {
    if (!file_.is_open()) {
        return;
    }

    std::string name = field_name;
    if (name.rfind("ELEMENT_NODAL_", 0) == 0) {
        name.erase(0, 14);
    }

    std::string base;
    const int frame = frame_number(name, base);
    const int key   = record_key(base);

    if (key == 0) return;

    const bool tensor = key == 11 || key == 21;
    const bool nodal  = field.domain == model::FieldDomain::NODE;

    if (!tensor && !nodal) return;
    if ( tensor && field.components != 6) return;
    if ( tensor && !nodal && field.domain != model::FieldDomain::ELEMENT_NODAL) return;

    logging::error(model_data_written_,
        "FilWriter: field '", field_name, "' written before the model definition");

    if (!increment_open_ || frame != current_frame_) {
        if (increment_open_) {
            end_increment();
        }
        begin_increment(frame, frame_value);
    }

    // Output request definition: 1 = nodal output, 0 = element output
    begin_record(1911);
    put_int (tensor ? 0 : 1);
    put_text("");
    put_text("");
    end_record();

    if (tensor) {
        logging::error(model_data != nullptr,
            "FilWriter: element field '", field_name, "' requires model data");
        write_element_field(field, key, *model_data);
    } else {
        write_nodal_field(field, key);
    }
}

/**
 * Writes records 1921 (release, date, sizes) and 1922 (heading).
 */
void FilWriter::write_header(const model::ModelData& model_data) {
    const auto& positions = *model_data.positions;

    const Index element_count = static_cast<Index>(
        std::count_if(element_types_.begin(), element_types_.end(),
                      [](const std::string& type) { return !type.empty(); }));

    // Typical element length estimated from the bounding box
    double diagonal = 0;
    for (Index d = 0; d < 3 && positions.rows > 0; ++d) {
        Precision lo = positions(0, d);
        Precision hi = positions(0, d);
        for (Index row = 1; row < positions.rows; ++row) {
            lo = std::min(lo, positions(row, d));
            hi = std::max(hi, positions(row, d));
        }
        diagonal += (hi - lo) * (hi - lo);
    }
    const double typical = std::sqrt(diagonal) / std::cbrt(static_cast<double>(std::max<Index>(element_count, 1)));

    char date[32] = {};
    char time[16] = {};
    const std::time_t now = std::time(nullptr);
    std::strftime(date, sizeof(date), "%d-%b-%Y", std::localtime(&now));
    std::strftime(time, sizeof(time), "%H:%M:%S", std::localtime(&now));

    begin_record(1921);
    put_text("FEMaster");
    put_text(date, 2);
    put_text(time);
    put_int (static_cast<std::int64_t>(element_count));
    put_int (static_cast<std::int64_t>(positions.rows));
    put_real(typical);
    end_record();

    begin_record(1922);
    put_text("generated by FEMaster", 10);
    end_record();
}

/**
 * Writes one 1901 record per compiled node, including isolated nodes.
 */
void FilWriter::write_nodes(const model::ModelData& model_data) {
    const auto& positions = *model_data.positions;

    for (Index row = 0; row < positions.rows; ++row) {
        begin_record(1901);
        put_int (node_ids_[row]);
        put_real(positions(row, 0));
        put_real(positions(row, 1));
        put_real(positions(row, 2));
        end_record();
    }
}

/**
 * Writes one 1900 record per element with an Abaqus equivalent. Every mapped
 * FEMaster element stores its nodes in the Abaqus input order of that type,
 * so connectivity is written unchanged.
 */
void FilWriter::write_elements(const model::ModelData& model_data) {
    for (std::size_t i = 0; i < model_data.elements.size(); ++i) {
        const auto& element = model_data.elements[i];
        if (!element || element_types_[i].empty()) continue;

        begin_record(1900);
        put_int (element_ids_[i]);
        put_text(element_types_[i]);
        for (ID node : *element) {
            put_int(node_ids_[static_cast<std::size_t>(node)]);
        }
        end_record();
    }
}

/**
 * Writes record 1902 from the complete structural DOF activation. A DOF is
 * active if any node activates it; active DOFs are numbered consecutively and
 * inactive DOFs are marked with zero. Purely thermal models fall back to the
 * three translations.
 */
void FilWriter::write_dofs(const SystemDofs& structural_dofs) {
    logging::error(structural_dofs.rows() == static_cast<Eigen::Index>(node_ids_.size())
                && structural_dofs.cols() == 6,
        "FilWriter: structural DOF mask does not match the compiled nodes");

    const bool any = structural_dofs.any();

    int count = 0;
    for (int d = 0; d < 6; ++d) {
        const bool active = any ? structural_dofs.col(d).any() : d < 3;
        active_dofs_[d]   = active ? ++count : 0;
    }

    begin_record(1902);
    for (int d = 0; d < 6; ++d) {
        put_int(active_dofs_[d]);
    }
    end_record();
}

/**
 * Opens an increment with record 2000. Eigenfrequency and buckling frames
 * additionally carry their eigenvalue in record 1980.
 *
 * Procedure keys: 1 static, 12 implicit dynamic, 41 frequency, 42 buckle.
 */
void FilWriter::begin_increment(int frame, Precision frame_value) {
    const bool   eigen        = current_step_type_ == WriterStepType::Eigenfrequency;
    const bool   buckle       = current_step_type_ == WriterStepType::Buckling;
    const double value        = std::isfinite(frame_value) ? frame_value : static_cast<double>(frame);
    const double time         = eigen || buckle ? 0.0 : value;
    const int    procedure    = eigen  ? 41
                              : buckle ? 42
                              : current_step_type_ == WriterStepType::Dynamic ? 12 : 1;

    begin_record(2000);
    put_real(time);
    put_real(time);
    put_real(0.0);
    put_real(0.0);
    put_int (procedure);
    put_int (current_step_);
    put_int (frame);
    put_int (eigen || buckle ? 1 : 0);
    put_real(0.0);
    put_real(eigen ? value : 0.0);
    put_real(0.0);
    put_text("FEMaster step " + std::to_string(current_step_), 10);
    end_record();

    // Eigenvalue record: omega^2 for frequency steps, load factor for buckling
    if (eigen || buckle) {
        const double omega = 2.0 * math::pi * value;

        begin_record(1980);
        put_int (frame);
        put_real(eigen ? omega * omega : value);
        if (eigen) {
            for (int i = 0; i < 2 + 12; ++i) put_real(0.0);
        }
        end_record();
    }

    increment_open_ = true;
    current_frame_  = frame;
}

/**
 * Closes the increment with record 2001, zero-padded so the block ends exactly.
 */
void FilWriter::end_increment() {
    Index remaining = kBlockWords - block_fill_;
    if (remaining < 2) {
        remaining += kBlockWords;
    }

    begin_record(2001);
    for (Index i = 2; i < remaining; ++i) put_int(0);
    end_record();

    increment_open_ = false;
}

/**
 * Writes one nodal record per node holding at least one finite component.
 * DOF fields are reduced to the active DOFs of record 1902; missing or
 * non-finite components are written as zero.
 */
void FilWriter::write_nodal_field(const model::Field& field,
                                  int                 key) {
    logging::error(node_ids_.size() == static_cast<std::size_t>(field.rows),
        "FilWriter: node id mapping does not match field rows");

    for (Index row = 0; row < field.rows; ++row) {
        bool any = false;
        for (Index c = 0; c < field.components; ++c) any = any || std::isfinite(field(row, c));
        if (!any) continue;

        begin_record(key);
        put_int(node_ids_[row]);
        if (key == 201) {
            put_real(std::isfinite(field(row, 0)) ? field(row, 0) : 0.0);
        } else {
            for (Index d = 0; d < 6; ++d) {
                if (active_dofs_[d] == 0) continue;
                const double v = d < field.components ? field(row, d) : 0.0;
                put_real(std::isfinite(v) ? v : 0.0);
            }
        }
        end_record();
    }
}

/**
 * Writes a tensor field as element output. NODE fields become values averaged
 * at nodes (location 4), ELEMENT_NODAL fields element-node values (location 2)
 * for every element present in the mesh.
 */
void FilWriter::write_element_field(const model::Field&     field,
                                    int                     key,
                                    const model::ModelData& model_data) {
    if (field.domain == model::FieldDomain::NODE) {
        for (Index row = 0; row < field.rows; ++row) {
            write_tensor_row(field, row, key, node_ids_[row], 0, 4);
        }
        return;
    }

    const auto& offsets = model_data.element_nodal_offsets;

    logging::error(offsets != nullptr,
        "FilWriter: element-nodal offsets are not initialized");

    for (std::size_t i = 0; i < model_data.elements.size(); ++i) {
        const auto& element = model_data.elements[i];
        if (!element || element_types_[i].empty()) continue;

        const Index offset = static_cast<Index>((*offsets)(i, 0));
        for (Index local = 0; local < element->n_nodes(); ++local) {
            const std::int64_t node = node_ids_[static_cast<std::size_t>(element->nodes()[local])];
            write_tensor_row(field, offset + local, key, element_ids_[i], node, 2);
        }
    }
}

/**
 * Writes the element header record 1 followed by one tensor data record in
 * Abaqus component order [11, 22, 33, 12, 13, 23].
 */
void FilWriter::write_tensor_row(const model::Field& field,
                                 Index               row,
                                 int                 key,
                                 std::int64_t        label,
                                 std::int64_t        point,
                                 int                 location) {
    if (!field.has_valid_row(row)) return;

    begin_record(1);
    put_int (label);
    put_int (point);
    put_int (0);
    put_int (location);
    put_text("");
    put_int (3);
    put_int (3);
    put_int (0);
    put_int (0);
    end_record();

    begin_record(key);
    for (Index c : {0, 1, 2, 5, 4, 3}) {
        put_real(field(row, c));
    }
    end_record();
}

/**
 * Maps a FEMaster element type name to its Abaqus element library name.
 * Returns an empty string for types without an Abaqus equivalent.
 */
std::string FilWriter::abaqus_type(const std::string& type_name) {
    if (type_name == "C3D4"   || type_name == "C3D6"  || type_name == "C3D8"
     || type_name == "C3D8R"  || type_name == "C3D8I" || type_name == "C3D10"
     || type_name == "C3D15"  || type_name == "C3D20" || type_name == "C3D20R"
     || type_name == "B31"    || type_name == "B33")  return type_name;
    if (type_name == "T3")                            return "T3D2";
    if (type_name == "MITC3FRT")                      return "S3";
    if (type_name == "MITC4FRT")                      return "S4";
    if (type_name == "MITC6FRT")                      return "STRI65";
    if (type_name == "MITC8FRT")                      return "S8R";
    return {};
}

/**
 * Maps a FEMaster base field name to its Abaqus record key, 0 if unsupported.
 */
int FilWriter::record_key(const std::string& base_name) {
    if (base_name == "DISPLACEMENT"
     || base_name == "MODE_SHAPE"
     || base_name == "BUCKLING_MODE")   return 101;
    if (base_name == "VELOCITY")        return 102;
    if (base_name == "ACCELERATION")    return 103;
    if (base_name == "REACTION_FORCES") return 104;
    if (base_name == "EXTERNAL_FORCES") return 106;
    if (base_name == "TEMPERATURE")     return 201;
    if (base_name == "STRESS")          return 11;
    if (base_name == "STRAIN")          return 21;
    return 0;
}

/**
 * Splits a trailing `_<n>` frame suffix from a field name. Fields without a
 * numeric suffix belong to frame one.
 */
int FilWriter::frame_number(const std::string& field_name, std::string& base_name) {
    const std::size_t split = field_name.rfind('_');
    const bool        tail  = split != std::string::npos && split + 1 < field_name.size()
                           && std::all_of(field_name.begin() + split + 1, field_name.end(),
                                          [](unsigned char c) { return std::isdigit(c); });

    base_name = tail ? field_name.substr(0, split) : field_name;
    return      tail ? std::stoi(field_name.substr(split + 1)) : 1;
}

/**
 * Starts a record with a length placeholder and its key.
 */
void FilWriter::begin_record(int key) {
    record_.clear();
    record_.push_back(0);
    put_int(key);
}

void FilWriter::put_int(std::int64_t value) {
    std::uint64_t word;
    std::memcpy(&word, &value, sizeof(word));
    record_.push_back(word);
}

void FilWriter::put_real(double value) {
    std::uint64_t word;
    std::memcpy(&word, &value, sizeof(word));
    record_.push_back(word);
}

/**
 * Appends text as `words` A8 words, left-justified and blank-padded.
 */
void FilWriter::put_text(const std::string& text, Index words) {
    std::string padded = text.substr(0, words * 8);
    padded.resize(words * 8, ' ');

    for (Index w = 0; w < words; ++w) {
        std::uint64_t word;
        std::memcpy(&word, padded.data() + w * 8, sizeof(word));
        record_.push_back(word);
    }
}

/**
 * Stores the record length in word one and moves the record into the stream.
 */
void FilWriter::end_record() {
    const std::int64_t length = static_cast<std::int64_t>(record_.size());
    std::memcpy(&record_[0], &length, sizeof(length));

    for (std::uint64_t word : record_) {
        put_word(word);
    }
    record_.clear();
}

void FilWriter::put_word(std::uint64_t word) {
    block_[block_fill_++] = word;
    if (block_fill_ == kBlockWords) {
        flush_block();
    }
}

/**
 * Writes the pending block as one Fortran unformatted record, zero-padded to
 * 512 words. Empty blocks are not written.
 */
void FilWriter::flush_block() {
    if (block_fill_ == 0) return;

    std::fill(block_.begin() + block_fill_, block_.end(), 0);

    const std::int32_t marker = static_cast<std::int32_t>(kBlockWords * sizeof(std::uint64_t));
    file_.write(reinterpret_cast<const char*>(&marker),       sizeof(marker));
    file_.write(reinterpret_cast<const char*>(block_.data()), static_cast<std::streamsize>(marker));
    file_.write(reinterpret_cast<const char*>(&marker),       sizeof(marker));

    block_fill_ = 0;
}

} // namespace writer
} // namespace io
} // namespace fem
