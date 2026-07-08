#pragma once

#include "writer_frd.h"
#include "writer_res.h"
#include "writer_file_formats.h"

#include <memory>
#include <string>

namespace fem {
namespace model {
struct ModelData;
struct Field;
}

namespace reader {

/**
 * @brief Facade for writing multiple result formats at once.
 *
 * Solver/loadcase code should only talk to this class.
 */
class ResultWriters {
    private:
    std::unique_ptr<ResWriter> res_writer;
    std::unique_ptr<FrdWriter> frd_writer;

    public:
    ResultWriters(const std::string& job_base_name,
                  const WriterFileFormats& options = WriterFileFormats());

    ResultWriters(ResultWriters&& other) noexcept;
    ResultWriters& operator=(ResultWriters&& other) noexcept;

    ResultWriters(const ResultWriters&) = delete;
    ResultWriters& operator=(const ResultWriters&) = delete;

    ~ResultWriters();

    void close();

    void write_model_data(const model::ModelData& model_data);

    void add_loadcase(int id);

    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr);
};

} // namespace reader
} // namespace fem