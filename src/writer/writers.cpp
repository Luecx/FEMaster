#include "writers.h"

#include "../data/field.h"
#include "../model/model_data.h"

#include <utility>

namespace fem {
namespace reader {

ResultWriters::ResultWriters(const std::string& job_base_name,
                             const WriterFileFormats& options) {
    if (job_base_name.empty()) {
        return;
    }

    if (options.res) {
        res_writer.reset(new ResWriter(job_base_name + ".res"));
    }

    if (options.frd) {
        frd_writer.reset(new FrdWriter(job_base_name + ".frd"));
    }
}

ResultWriters::~ResultWriters() {
    close();
}

ResultWriters::ResultWriters(ResultWriters&& other) noexcept
    : res_writer(std::move(other.res_writer)),
      frd_writer(std::move(other.frd_writer)) {}

ResultWriters& ResultWriters::operator=(ResultWriters&& other) noexcept {
    if (this != &other) {
        close();

        res_writer = std::move(other.res_writer);
        frd_writer = std::move(other.frd_writer);
    }

    return *this;
}

void ResultWriters::close() {
    if (res_writer) {
        res_writer->close();
    }

    if (frd_writer) {
        frd_writer->close();
    }
}

void ResultWriters::write_model_data(const model::ModelData& model_data) {
    if (frd_writer) {
        frd_writer->write_model_data(model_data);
    }
}

void ResultWriters::add_loadcase(int id) {
    if (res_writer) {
        res_writer->add_loadcase(id);
    }

    if (frd_writer) {
        frd_writer->add_loadcase(id);
    }
}

void ResultWriters::write_field(const model::Field& field,
                                const std::string& field_name,
                                const model::ModelData* model_data) {
    if (res_writer) {
        res_writer->write_field(field, field_name, model_data);
    }

    if (frd_writer) {
        frd_writer->write_field(field, field_name, model_data);
    }
}

} // namespace reader
} // namespace fem
