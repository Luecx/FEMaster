/**
 * @file writer_fil.h
 * @brief Declares the binary Abaqus results-file (.fil) writer.
 *
 * `FilWriter` emits the compiled FEMaster mesh and supported result fields in
 * the sequential Abaqus/Standard results-file format. Every record consists of
 * 8-byte words
 *
 *     [length] [key] [attribute 1] ... [attribute n],
 *
 * where integers, reals and A8 character strings each occupy one word. Words
 * form a continuous stream that is cut into blocks of 512 words, each written
 * as one Fortran unformatted sequential record (4-byte length markers before
 * and after the 4096 data bytes). Records may cross block boundaries.
 *
 * Every increment starts with record 2000 and ends with record 2001. Record
 * 2001 is padded with zero words so that the next increment starts on a new
 * block, as required by Abaqus' FILE utilities.
 *
 * Node and element labels follow the FRD convention
 *
 *     10^8 * i_instance + n_local,
 *
 * so FRD and FIL output of the same job reference identical entities.
 *
 * Unsupported fields (step summaries, topology optimization fields, ...) are
 * ignored intentionally.
 *
 * @see FrdWriter
 * @see model::ModelData
 *
 * @author Finn Eggers
 * @date 17.11.2026
 */

#pragma once

#include "../../core/core.h"
#include "../../core/logging.h"
#include "../../data/field.h"
#include "writer_step_type.h"

#include <array>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
}

namespace io {
namespace writer {

/**
 * @brief Binary Abaqus `.fil` mesh, nodal and element result writer.
 *
 * Supported result records:
 *
 *  | FEMaster field                          | key | Abaqus |
 *  |-----------------------------------------|-----|--------|
 *  | DISPLACEMENT, MODE_SHAPE, BUCKLING_MODE | 101 | U      |
 *  | VELOCITY                                | 102 | V      |
 *  | ACCELERATION                            | 103 | A      |
 *  | REACTION_FORCES                         | 104 | RF     |
 *  | EXTERNAL_FORCES                         | 106 | CF     |
 *  | TEMPERATURE                             | 201 | NT     |
 *  | STRESS                                  |  11 | S      |
 *  | STRAIN                                  |  21 | E      |
 *
 * Nodal records contain the active degrees of freedom announced in record
 * 1902. Tensor fields are written as element output (record 1 header followed
 * by the data record) at NODE (location 4, averaged at nodes), ELEMENT_NODAL
 * (location 2) or ELEMENT_IP (location 0). FEMaster's Voigt order
 * [xx, yy, zz, yz, zx, xy] is permuted to Abaqus' [11, 22, 33, 12, 13, 23];
 * both use engineering shear strains.
 */
class FilWriter {
    private:
    static constexpr Index kBlockWords = 512;

    std::ofstream file_;

    // Dense ids -> external labels, active-DOF positions (0 = inactive)
    std::vector<std::int64_t> node_ids_;
    std::vector<std::int64_t> element_ids_;
    std::array<int, 6>        active_dofs_{};
    int                       active_count_ = 0;

    // Word stream: record under construction and the pending 512-word block
    std::vector<std::uint64_t>             record_;
    std::array<std::uint64_t, kBlockWords> block_{};
    Index                                  block_fill_ = 0;

    // Writer state
    bool           model_data_written_ = false;
    bool           increment_open_     = false;
    int            current_step_       = 1;
    int            current_frame_      = -1;
    WriterStepType current_step_type_  = WriterStepType::Static;

    public:
    // File lifetime and ownership
    explicit FilWriter(const std::string& filename = "");
    ~FilWriter();

    FilWriter(const FilWriter&)            = delete;
    FilWriter& operator=(const FilWriter&) = delete;

    void open (const std::string& filename);
    void close();

    // Analysis-step metadata
    void add_loadcase(int id, WriterStepType step_type = WriterStepType::Static);

    // Mesh and result output
    void write_model_data(const model::ModelData& model_data);
    void write_field     (const model::Field&     field,
                          const std::string&      field_name,
                          const model::ModelData* model_data  = nullptr,
                          Precision               frame_value = std::numeric_limits<Precision>::quiet_NaN());

    private:
    // Model definition records 1921, 1922, 1901, 1900, 1902
    void write_header  (const model::ModelData& model_data);
    void write_nodes   (const model::ModelData& model_data);
    void write_elements(const model::ModelData& model_data);
    void write_dofs    (const model::ModelData& model_data);

    // Increment records 2000, 1980, 2001
    void begin_increment(int frame, Precision frame_value);
    void end_increment  ();

    // Result records
    void write_nodal_field  (const model::Field&     field,
                             int                     key);
    void write_element_field(const model::Field&     field,
                             int                     key,
                             const model::ModelData& model_data);
    void write_tensor_row   (const model::Field&     field,
                             Index                   row,
                             int                     key,
                             std::int64_t            label,
                             std::int64_t            point,
                             int                     location);

    // FEMaster field name -> record key (0 if unsupported) and frame number
    static int record_key  (const std::string& base_name);
    static int frame_number(const std::string& field_name, std::string& base_name);

    // Record assembly and blocked word output
    void begin_record(int key);
    void put_int     (std::int64_t value);
    void put_real    (double value);
    void put_text    (const std::string& text, Index words = 1);
    void end_record  ();
    void put_word    (std::uint64_t word);
    void flush_block ();
};

} // namespace writer
} // namespace io
} // namespace fem
