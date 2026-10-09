/**
 * @file writer_fil.h
 * @brief Abaqus ASCII .fil writer (80-column sequential records).
 */
#pragma once

#include "../../core/core.h"
#include "../../data/field.h"
#include "writer_step_type.h"

#include <fstream>
#include <limits>
#include <string>
#include <string_view>

namespace fem::model { struct ModelData; }

namespace fem::io::writer {

class FilWriter {
public:
    explicit FilWriter(const std::string& filename = "");
    ~FilWriter();

    FilWriter(const FilWriter&) = delete;
    FilWriter& operator=(const FilWriter&) = delete;

    void open(const std::string& filename);
    void close();
    void write_model_data(const model::ModelData& model_data);
    void add_loadcase(int id, WriterStepType type = WriterStepType::Static);
    void begin_frame(Precision value = std::numeric_limits<Precision>::quiet_NaN());
    void end_frame();
    void write_field(const model::Field& field, const std::string& name,
                     const model::ModelData* data = nullptr,
                     Precision value = std::numeric_limits<Precision>::quiet_NaN());

private:
    void append(std::string_view value);
    void finish_line();
    void record(int key, std::size_t attributes);
    void integer(long long value);
    void real(double value);
    void string(std::string_view value);
    void write_nodal(const model::Field& field, int key);
    void write_element(const model::Field& field, int key);
    void element_header(int element, int point, int location);
    void output_request(bool nodal, std::string_view type = {});

    std::ofstream file_;
    std::string line_;
    const model::ModelData* model_data_ = nullptr;
    bool model_written_ = false;
    bool frame_open_ = false;
    bool rotations_ = false;
    int step_ = 1;
    int increment_ = 0;
    WriterStepType step_type_ = WriterStepType::Static;
    double previous_step_time_ = 0.;
    double total_time_offset_ = 0.;
};
} // namespace fem::io::writer
