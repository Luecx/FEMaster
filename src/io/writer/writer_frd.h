/**
 * @file writer_frd.h
 * @brief Declares the ASCII CalculiX/CGX FRD result writer.
 *
 * `FrdWriter` writes the compiled FEMaster mesh and NODE-domain result fields
 * in the CalculiX FRD format. Dense solver node identifiers are translated once
 * during mesh output into persistent FRD node identifiers derived from the
 * semantic Instance and part-local node identifiers retained by `ModelData`.
 *
 * The resulting dense-to-FRD node mapping is cached by the writer and reused
 * for element connectivity and all later nodal result blocks. Result output
 * therefore does not require repeated semantic node resolution.
 *
 * Non-NODE fields are ignored intentionally. Element identifiers are written
 * unchanged from the compiled model.
 *
 * @see FrdWriter
 * @see model::ModelData
 *
 * @author Finn Eggers
 * @date 20.08.2026
 */

#pragma once

#include "../../core/core.h"
#include "../../core/logging.h"
#include "../../data/field.h"
#include "writer_step_type.h"

#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
struct ElementInterface;
}

namespace io {
namespace writer {

/**
 * @brief Minimal ASCII CalculiX/CGX FRD mesh and nodal-result writer.
 *
 * The writer emits all compiled nodes, supported structural element
 * connectivities and NODE-domain result fields. Nodes are not filtered by
 * element participation, so isolated model nodes remain present in the output.
 *
 * FEMaster uses dense node identifiers internally. During the first mesh write,
 * each dense node is assigned its FRD identifier according to
 *
 *     10^8 * i_instance + n_local.
 *
 * The implicit default Instance has instance identifier zero, so its original
 * local node identifiers remain unchanged in the FRD file. Explicit Instances
 * occupy disjoint identifier ranges of width \f$10^8\f$.
 *
 * The generated dense-to-FRD mapping is retained by the writer for subsequent
 * field output, allowing result frames to be written without repeated access to
 * the complete model data.
 */
class FrdWriter {
    private:
    std::ofstream file_path;

    // Dense global node id -> FRD node id
    std::vector<ID> node_ids;

    // Writer state
    bool model_data_written = false;

    int current_step         = 1;
    int current_result_block = 0;
    int current_global_frame = 0;
    int last_frame_step      = -1;
    int last_frame_id        = -1;

    WriterStepType current_step_type = WriterStepType::Static;

    public:
    // File lifetime and ownership
    explicit FrdWriter(const std::string& filename = "");
    ~FrdWriter();

    FrdWriter(const FrdWriter&) = delete;
    FrdWriter& operator=(const FrdWriter&) = delete;

    void open(const std::string& filename);
    void close();

    // Analysis-step metadata
    void add_loadcase(int id, WriterStepType step_type = WriterStepType::Static);

    // Mesh and result output
    void write_model_data(const model::ModelData& model_data);
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr,
                     Precision frame_value = std::numeric_limits<Precision>::quiet_NaN());

    private:
    // FRD mesh and NODE-domain result blocks
    void write_nodes      (const model::ModelData& model_data);
    void write_elements   (const model::ModelData& model_data);
    void write_nodal_field(const model::Field& field,
                           const std::string& field_name,
                           Precision frame_value);
};

} // namespace writer
} // namespace io
} // namespace fem