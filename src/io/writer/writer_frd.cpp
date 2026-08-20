/**
 * @file writer_frd.cpp
 * @brief Implements the ASCII CalculiX/CGX FRD result writer.
 *
 * The writer emits compiled FEMaster nodes, supported structural elements and
 * NODE-domain result fields in CalculiX FRD syntax.
 *
 * FEMaster solver topology uses dense node identifiers. When model data is
 * written, the semantic reverse mapping retained by `ModelData` is evaluated
 * once to construct the external FRD identifier
 *
 * \f[
 *     n_\mathrm{FRD}
 *     =
 *     10^8\,i_\mathrm{instance}
 *     +
 *     n_\mathrm{local}.
 * \f]
 *
 * The resulting dense-to-FRD mapping is cached locally in `FrdWriter`. Node
 * coordinates, element connectivity and every subsequent nodal result block
 * therefore use identical external identifiers without repeated Instance
 * lookups.
 *
 * All compiled nodes are emitted, including nodes not referenced by supported
 * elements. Compiled element identifiers are written unchanged.
 *
 * @see FrdWriter
 * @see model::ModelData
 *
 * @author Finn Eggers
 * @date 20.08.2026
 */

#include "writer_frd.h"

#include "../../model/element/element.h"
#include "../../model/element/element_structural.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <vector>

namespace fem {
namespace io {
namespace writer {

namespace {

// One FRD result component, e.g. the stress Sxy component
struct FRDComponent {
    std::string name;
    int entity   = 1; // scalar = 1, vector = 2, tensor = 4
    int index_1  = 0; // first vector or tensor index
    int index_2  = 0; // second tensor index
    bool derived = false;

    static FRDComponent scalar(const std::string& name, int index   ) { return {name, 1, index, 0, false};}
    static FRDComponent vector(const std::string& name, int index   ) { return {name, 2, index, 0, false};}
    static FRDComponent tensor(const std::string& name, int i, int j) { return {name, 4, i, j, false};}
    static FRDComponent vcnorm(const std::string& name)               { return {name, 2, 0, 0, true };}
};

// FRD field definition and accepted FEMaster aliases
struct FRDField {
    std::vector<std::string> aliases;
    std::string name;
    std::vector<FRDComponent> components;
};

/**
 * Extracts the numeric frame identifier encoded in a field name.
 *
 * All decimal digits occurring in the field name are concatenated. Fields
 * without an explicit numeric suffix belong to frame one.
 *
 * @param s Field name supplied by the result writer.
 * @return Positive FRD frame identifier.
 */
int frd_frame(const std::string& s) {
    std::string out;

    for (unsigned char c : s) {
        if (std::isdigit(c)) {
            out += static_cast<char>(c);
        }
    }

    return out.empty() ? 1 : std::stoi(out);
}

/**
 * Maps one structural FEMaster element to its CalculiX FRD topology code.
 *
 * Shell and solid topology is identified from the structural-element category
 * together with its node count. Remaining two- and three-node structural
 * elements are interpreted as line elements. Unsupported or non-structural
 * elements return zero and are omitted from the FRD element block.
 *
 * @param element Compiled element to classify.
 * @return FRD topology code or zero if unsupported.
 */
int frd_element_type(const model::ElementInterface& element) {
    const auto* structural = dynamic_cast<const model::StructuralElement*>(&element);

    if (!structural) {
        return 0;
    }

    const Dim n_nodes = element.n_nodes();

    if (structural->is_shell()) {
        if (n_nodes == 3) return 7;
        if (n_nodes == 6) return 8;
        if (n_nodes == 4) return 9;
        if (n_nodes == 8) return 10;

        return 0;
    }

    if (structural->is_solid()) {
        if (n_nodes == 4)  return 3;
        if (n_nodes == 10) return 6;
        if (n_nodes == 8)  return 1;
        if (n_nodes == 20) return 4;
        if (n_nodes == 6)  return 2;
        if (n_nodes == 15) return 5;

        return 0;
    }

    if (n_nodes == 2) return 11;
    if (n_nodes == 3) return 12;

    return 0;
}

/**
 * Maps the FEMaster analysis-step category to the FRD increment type.
 *
 * @param step_type Analysis category of the current writer step.
 * @return CalculiX FRD increment-type identifier.
 */
int frd_increment_type(WriterStepType step_type) {
    switch (step_type) {
        case WriterStepType::Dynamic:        return 1;
        case WriterStepType::Eigenfrequency: return 2;
        case WriterStepType::Buckling:       return 4;
        case WriterStepType::Static:         return 0;
    }

    return 0;
}

/**
 * Returns the FRD analysis label associated with one FEMaster step type.
 *
 * Static results do not require an explicit analysis label and therefore return
 * an empty string.
 *
 * @param step_type Analysis category of the current writer step.
 * @return FRD analysis label.
 */
const char* frd_analysis_name(WriterStepType step_type) {
    switch (step_type) {
        case WriterStepType::Dynamic:        return "DYNAMIC";
        case WriterStepType::Eigenfrequency: return "MODAL";
        case WriterStepType::Buckling:       return "BUCKLING";
        case WriterStepType::Static:         return "";
    }

    return "";
}

/**
 * Resolves a FEMaster field name to its FRD field and component definition.
 *
 * Field aliases are matched case-insensitively after removing all nonalphabetic
 * characters. Component metadata describes scalar, vector and tensor semantics
 * required by the FRD result header.
 *
 * @param field_name FEMaster result-field name.
 * @return Matching immutable FRD field definition.
 */
const FRDField& frd_field(const std::string& field_name) {
    std::string name;

    for (unsigned char c : field_name) {
        if (std::isalpha(c)) {
            name += static_cast<char>(std::toupper(c));
        }
    }

    static const std::vector<FRDField> definitions{
        {
            {"DISPLACEMENT", "DISP", "MODESHAPE", "BUCKLINGMODE"}, "DISP",
            {
                FRDComponent::vector("D1" , 1),
                FRDComponent::vector("D2" , 2),
                FRDComponent::vector("D3" , 3),
                FRDComponent::scalar("D4" , 4),
                FRDComponent::scalar("D5" , 5),
                FRDComponent::scalar("D6" , 6),
                FRDComponent::vcnorm("ALL")
            }
        },
        {
            {"VELOCITY", "VELO"}, "VELO",
            {
                FRDComponent::vector("V1" , 1),
                FRDComponent::vector("V2" , 2),
                FRDComponent::vector("V3" , 3),
                FRDComponent::vcnorm("ALL")
            }
        },
        {
            {"ACCELERATION", "ACCE"}, "ACCE",
            {
                FRDComponent::vector("A1" , 1),
                FRDComponent::vector("A2" , 2),
                FRDComponent::vector("A3" , 3),
                FRDComponent::vcnorm("ALL")
            }
        },
        {
            {"REACTIONFORCES", "FORC"}, "FORC",
            {
                FRDComponent::vector("F1" , 1),
                FRDComponent::vector("F2" , 2),
                FRDComponent::vector("F3" , 3),
                FRDComponent::scalar("F4" , 4),
                FRDComponent::scalar("F5" , 5),
                FRDComponent::scalar("F6" , 6),
                FRDComponent::vcnorm("ALL")
            }
        },
        {
            {"EXTERNALFORCES", "EXTFORC"}, "EXTFORC",
            {
                FRDComponent::vector("F1" , 1),
                FRDComponent::vector("F2" , 2),
                FRDComponent::vector("F3" , 3),
                FRDComponent::scalar("F4" , 4),
                FRDComponent::scalar("F5" , 5),
                FRDComponent::scalar("F6" , 6),
                FRDComponent::vcnorm("ALL")
            }
        },
        {
            {"INTERNALFORCES", "INTFORC"}, "INTFORC",
            {
                FRDComponent::vector("F1" , 1),
                FRDComponent::vector("F2" , 2),
                FRDComponent::vector("F3" , 3),
                FRDComponent::scalar("F4" , 4),
                FRDComponent::scalar("F5" , 5),
                FRDComponent::scalar("F6" , 6),
                FRDComponent::vcnorm("ALL")
            }
        },
        {
            {"STRESS"}, "STRESS",
            {
                FRDComponent::tensor("SXX", 1, 1),
                FRDComponent::tensor("SYY", 2, 2),
                FRDComponent::tensor("SZZ", 3, 3),
                FRDComponent::tensor("SYZ", 2, 3),
                FRDComponent::tensor("SZX", 3, 1),
                FRDComponent::tensor("SXY", 1, 2)
            }
        },
        {
            {"STRAIN", "TOTALSTRAIN", "TOSTRAIN"}, "TOSTRAIN",
            {
                FRDComponent::tensor("EXX", 1, 1),
                FRDComponent::tensor("EYY", 2, 2),
                FRDComponent::tensor("EZZ", 3, 3),
                FRDComponent::tensor("EYZ", 2, 3),
                FRDComponent::tensor("EZX", 3, 1),
                FRDComponent::tensor("EXY", 1, 2)
            }
        },
        {
            {"MECHANICALSTRAIN", "MESTRAIN"}, "MESTRAIN",
            {
                FRDComponent::tensor("MEXX", 1, 1),
                FRDComponent::tensor("MEYY", 2, 2),
                FRDComponent::tensor("MEZZ", 3, 3),
                FRDComponent::tensor("MEYZ", 2, 3),
                FRDComponent::tensor("MEZX", 3, 1),
                FRDComponent::tensor("MEXY", 1, 2)
            }
        },
        {
            {"STRESSTOP"}, "STRPOS",
            {
                FRDComponent::tensor("SXX", 1, 1),
                FRDComponent::tensor("SYY", 2, 2),
                FRDComponent::tensor("SZZ", 3, 3),
                FRDComponent::tensor("SYZ", 2, 3),
                FRDComponent::tensor("SZX", 3, 1),
                FRDComponent::tensor("SXY", 1, 2)
            }
        },
        {
            {"STRESS"}, "STRMID",
            {
                FRDComponent::tensor("SXX", 1, 1),
                FRDComponent::tensor("SYY", 2, 2),
                FRDComponent::tensor("SZZ", 3, 3),
                FRDComponent::tensor("SYZ", 2, 3),
                FRDComponent::tensor("SZX", 3, 1),
                FRDComponent::tensor("SXY", 1, 2)
            }
        },
        {
            {"STRESSBOT"}, "STRNEG",
            {
                FRDComponent::tensor("SXX", 1, 1),
                FRDComponent::tensor("SYY", 2, 2),
                FRDComponent::tensor("SZZ", 3, 3),
                FRDComponent::tensor("SYZ", 2, 3),
                FRDComponent::tensor("SZX", 3, 1),
                FRDComponent::tensor("SXY", 1, 2)
            }
        },
        {
            {"SHELLRESULTANTS"}, "SHR",
            {
                FRDComponent::vector("NXX", 1),
                FRDComponent::vector("NYY", 2),
                FRDComponent::vector("NXY", 3),
                FRDComponent::vector("MXX", 4),
                FRDComponent::vector("MYY", 5),
                FRDComponent::vector("MXY", 6),
                FRDComponent::vector("QX", 7),
                FRDComponent::vector("QY", 8)
            }
        },
    };

    for (const FRDField& definition : definitions) {
        if (std::find(definition.aliases.begin(),
                      definition.aliases.end(),
                      name) != definition.aliases.end()) {
            return definition;
        }
    }

    logging::error(false,
        "FrdWriter: unsupported NODE field name: ", field_name);

    return definitions.front();
}

/**
 * Returns one finite FRD field component value.
 *
 * Non-finite stored values are written as zero. Derived vector norms use the
 * first three field components and replace non-finite individual components by
 * zero before evaluating their Euclidean norm.
 *
 * @param field Source FEMaster field.
 * @param row Nodal field row.
 * @param component Component index for directly stored values.
 * @param frd_component FRD component metadata.
 * @return Finite scalar value for output.
 */
Precision field_value(const model::Field& field,
                      Index row,
                      Index component,
                      const FRDComponent& frd_component) {
    if (frd_component.derived) {
        const Precision x = std::isfinite(field(row, 0)) ? field(row, 0) : 0.0;
        const Precision y = std::isfinite(field(row, 1)) ? field(row, 1) : 0.0;
        const Precision z = std::isfinite(field(row, 2)) ? field(row, 2) : 0.0;

        return std::sqrt(x * x + y * y + z * z);
    }

    return std::isfinite(field(row, component))
        ? field(row, component)
        : 0;
}

/**
 * Writes one finite floating-point value in the fixed FRD scientific format.
 *
 * Non-finite values are replaced by zero to prevent invalid numerical tokens
 * from entering the result file.
 *
 * @param file_path Active FRD output stream.
 * @param value Value to write.
 * @param width Output field width.
 */
void write_float(std::ofstream& file_path,
                 Precision value,
                 int width = 12) {
    file_path << std::right
              << std::uppercase
              << std::scientific
              << std::setprecision(5)
              << std::setw(width)
              << (std::isfinite(value) ? value : 0);
}

} // namespace

/**
 * Constructs an FRD writer and optionally opens its output file immediately.
 *
 * @param filename Output path or an empty string to construct a closed writer.
 */
FrdWriter::FrdWriter(const std::string& filename) {
    if (!filename.empty()) {
        open(filename);
    }
}

/**
 * Finalizes and closes the active FRD file.
 */
FrdWriter::~FrdWriter() {
    close();
}

/**
 * Opens a new FRD file and resets all writer-local state.
 *
 * Any previously open file is finalized first. The dense-to-FRD node mapping is
 * cleared because it belongs to the mesh written into one specific output file.
 *
 * @param filename Output path to open with truncation.
 */
void FrdWriter::open(const std::string& filename) {
    close();

    file_path.open(filename, std::ios::out | std::ios::trunc);

    logging::error(file_path.is_open(),
        "FrdWriter: failed to open file: ", filename);

    // Reset mesh numbering and analysis-frame state
    node_ids.clear();

    model_data_written   = false;
    current_step         = 1;
    current_result_block = 0;
    current_global_frame = 0;
    last_frame_step      = -1;
    last_frame_id        = -1;
    current_step_type    = WriterStepType::Static;

    file_path << "    1CFEMAST\n";
    file_path << "    1Ugenerated by FEMaster\n";
}

/**
 * Finalizes the active FRD stream and closes the file.
 *
 * A closed writer is left untouched.
 */
void FrdWriter::close() {
    if (file_path.is_open()) {
        file_path << " 9999\n";
        file_path.close();
    }
}

/**
 * Selects the analysis-step metadata used by subsequent result blocks.
 *
 * Non-positive loadcase identifiers fall back to one. Frame tracking is reset
 * while the global result-frame counter remains monotonic across loadcases.
 *
 * @param id Analysis-step identifier.
 * @param step_type Physical type of the analysis step.
 */
void FrdWriter::add_loadcase(int id, WriterStepType step_type) {
    logging::error(file_path.is_open(),
        "FrdWriter: cannot add loadcase: file is not open");

    current_step      = id > 0 ? id : 1;
    current_step_type = step_type;
    last_frame_step   = -1;
    last_frame_id     = -1;
}

/**
 * Writes the compiled mesh and establishes its persistent FRD node numbering.
 *
 * Every dense FEMaster node is resolved once through `ModelData::node_mapping`.
 * The FRD identifier is formed from the owning Instance and the original
 * part-local node identifier as
 *
 * \f[
 *     n_\mathrm{FRD}
 *     =
 *     10^8\,i_\mathrm{instance}
 *     +
 *     n_\mathrm{local}.
 * \f]
 *
 * The local identifier must remain below \f$10^8\f$ so adjacent Instance ranges
 * cannot overlap. The completed numbering is cached and subsequently reused by
 * coordinate, connectivity and nodal-result output.
 *
 * @param model_data Compiled model containing topology, positions and semantic
 * node mappings.
 */
void FrdWriter::write_model_data(const model::ModelData& model_data) {
    logging::error(file_path.is_open(),
        "FrdWriter: cannot write model data: file is not open");
    logging::error(model_data.positions != nullptr,
        "FrdWriter: model_data.positions is not initialized");
    logging::error(model_data.node_mapping.size()
                    == static_cast<std::size_t>(model_data.positions->rows),
        "FrdWriter: node mapping does not match positions field");

    if (model_data_written) {
        return;
    }

    // Build the dense-to-FRD node numbering once for this output file
    node_ids.resize(model_data.node_mapping.size());

    for (std::size_t i = 0; i < model_data.node_mapping.size(); ++i) {
        const auto& [instance, local_id] = model_data.node_mapping[i];

        logging::error(instance != nullptr,
            "FrdWriter: dense node ", i, " has no source instance");
        logging::error(local_id < static_cast<ID>(100000000),
            "FrdWriter: local node id ", local_id, " exceeds the supported instance range");

        node_ids[i] =
              static_cast<Index>(100000000)
            * static_cast<Index>(instance->instance_id)
            + static_cast<Index>(local_id);
    }

    // Emit geometry and connectivity using the fixed external numbering
    write_nodes(model_data);
    write_elements(model_data);

    model_data_written = true;
}

/**
 * Writes one supported NODE-domain field.
 *
 * Model data is needed only when the mesh has not yet been written. Once the
 * mesh block has established the cached FRD node numbering, all later result
 * frames can be emitted from the field alone.
 *
 * Non-NODE fields and writes to a closed stream are ignored intentionally.
 *
 * @param field Result field to write.
 * @param field_name Semantic field name used to select its FRD representation.
 * @param model_data Model data required only for the initial mesh write.
 * @param frame_value Physical frame value or NaN to derive it from the field name.
 */
void FrdWriter::write_field(const model::Field& field,
                            const std::string& field_name,
                            const model::ModelData* model_data,
                            Precision frame_value) {
    if (!file_path.is_open()) {
        return;
    }

    if (field.domain != model::FieldDomain::NODE) {
        return;
    }

    if (!model_data_written) {
        logging::error(model_data != nullptr,
            "FrdWriter: NODE field '", field_name,
            "' requires model data because mesh was not written yet");

        write_model_data(*model_data);
    }

    write_nodal_field(field, field_name, frame_value);
}

/**
 * Writes the complete compiled nodal coordinate block.
 *
 * Every row of the compiled position field is emitted. Nodes are deliberately
 * not filtered by element connectivity, so isolated and otherwise unused nodes
 * remain represented in the FRD mesh.
 *
 * @param model_data Compiled model providing the nodal position field.
 */
void FrdWriter::write_nodes(const model::ModelData& model_data) {
    const auto& positions = *model_data.positions;

    logging::error(node_ids.size() == static_cast<std::size_t>(positions.rows),
        "FrdWriter: node id mapping does not match positions field");

    file_path << "    2C"
              << std::string(18, ' ')
              << std::setw(12)
              << static_cast<long long>(positions.rows)
              << std::string(37, ' ')
              << 1
              << '\n';

    for (Index row = 0; row < positions.rows; ++row) {
        file_path << " -1"
                  << std::setw(10)
                  << static_cast<long long>(
                         node_ids[static_cast<std::size_t>(row)]);

        write_float(file_path, positions(row, 0));
        write_float(file_path, positions(row, 1));
        write_float(file_path, positions(row, 2));

        file_path << '\n';
    }

    file_path << " -3\n";
}

/**
 * Writes all supported compiled structural elements and their connectivity.
 *
 * Element identifiers are emitted unchanged. Their connectivity contains dense
 * FEMaster node identifiers internally and is translated through the cached
 * dense-to-FRD node mapping immediately before output.
 *
 * CalculiX requires different connectivity ordering for 20-node hexahedra and
 * 15-node wedges, so these two element families are reordered explicitly.
 *
 * @param model_data Compiled model containing structural element topology.
 */
void FrdWriter::write_elements(const model::ModelData& model_data) {
    std::size_t element_count = 0;

    // Determine the number advertised in the FRD element block header
    for (const auto& element : model_data.elements) {
        if (element && frd_element_type(*element) != 0) {
            ++element_count;
        }
    }

    logging::error(element_count > 0,
        "FrdWriter: no supported elements found for FRD output");

    file_path << "    3C"
              << std::string(18, ' ')
              << std::setw(12)
              << static_cast<long long>(element_count)
              << std::string(37, ' ')
              << 1
              << '\n';

    // Emit every supported element in compiled order
    for (const auto& element : model_data.elements) {
        if (!element) {
            continue;
        }

        const int type = frd_element_type(*element);

        if (type == 0) {
            continue;
        }

        file_path << " -1"
                  << std::setw(10)
                  << static_cast<long long>(element->elem_id)
                  << std::setw(5) << type
                  << std::setw(5) << 0
                  << std::setw(5) << 0
                  << '\n';

        // Copy connectivity so FRD-specific high-order permutations remain local
        std::vector<ID> connectivity;
        connectivity.reserve(static_cast<std::size_t>(element->n_nodes()));

        for (ID node_id : *element) {
            connectivity.push_back(node_id);
        }

        if (type == 4 && connectivity.size() == 20) {
            const std::vector<ID> internal = connectivity;
            const int order[] = {
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                10, 11, 16, 17, 18, 19, 12, 13, 14, 15
            };

            for (std::size_t i = 0; i < connectivity.size(); ++i) {
                connectivity[i] = internal[order[i]];
            }
        }

        if (type == 5 && connectivity.size() == 15) {
            const std::vector<ID> internal = connectivity;
            const int order[] = {
                0, 1, 2, 3, 4, 5, 6, 7, 8, 12,
                13, 14, 9, 10, 11
            };

            for (std::size_t i = 0; i < connectivity.size(); ++i) {
                connectivity[i] = internal[order[i]];
            }
        }

        // Translate dense connectivity to the external FRD node numbering
        Dim local = 0;

        for (ID node_id : connectivity) {
            if (local % static_cast<Dim>(10) == 0) {
                if (local != 0) {
                    file_path << '\n';
                }

                file_path << " -2";
            }

            logging::error(static_cast<std::size_t>(node_id) < node_ids.size(),
                "FrdWriter: element ", element->elem_id,
                " references node ", node_id,
                " outside the FRD node mapping");

            file_path << std::setw(10)
                      << static_cast<long long>(
                             node_ids[static_cast<std::size_t>(node_id)]);

            ++local;
        }

        file_path << '\n';
    }

    file_path << " -3\n";
}

/**
 * Writes one NODE-domain result block using the cached FRD node numbering.
 *
 * Field rows correspond directly to dense global node identifiers. Each row is
 * therefore translated by indexing the same `node_ids` vector used for mesh
 * connectivity, guaranteeing consistent numbering across geometry and results.
 *
 * Derived vector norms are evaluated during scalar output through `field_value`.
 *
 * @param field Nodal result field.
 * @param field_name Semantic name selecting the FRD field definition.
 * @param frame_value Physical frame value or NaN to derive it from the field name.
 */
void FrdWriter::write_nodal_field(const model::Field& field,
                                  const std::string& field_name,
                                  Precision frame_value) {
    logging::error(node_ids.size() == static_cast<std::size_t>(field.rows),
        "FrdWriter: node id mapping does not match field rows");

    // Resolve field metadata and analysis-frame numbering
    const FRDField& frd = frd_field(field_name);
    int frame           = frd_frame(field_name);
    const int step      = current_step > 0 ? current_step : 1;
    const int block     = ++current_result_block;

    if (current_step_type == WriterStepType::Dynamic) {
        ++frame;
    }

    if (frame < 1) {
        frame = 1;
    }

    const Precision value = std::isfinite(frame_value)
        ? frame_value
        : static_cast<Precision>(frame);

    if (last_frame_step != step || last_frame_id != frame) {
        ++current_global_frame;
        last_frame_step = step;
        last_frame_id   = frame;
    }

    const int global_frame = current_global_frame;
    const int ictype       = frd_increment_type(current_step_type);
    const char* analysis   = frd_analysis_name(current_step_type);

    // Write FRD result-step metadata
    file_path << "    1PSTEP"
              << std::setw(26) << block
              << std::setw(12) << frame
              << std::setw(12) << step
              << '\n';

    if (current_step_type == WriterStepType::Eigenfrequency) {
        file_path << "    1PMODE"
                  << std::setw(26) << frame
                  << '\n';
    }

    file_path << "  100CL  "
              << std::setw(3) << 100 + global_frame
              << ' ';

    write_float(file_path, value, 11);

    file_path << std::setw(12) << field.rows
              << std::string(21, ' ')
              << ictype
              << std::setw(5) << global_frame;

    if (analysis[0] == '\0') {
        file_path << std::string(11, ' ');
    } else {
        file_path << std::left
                  << std::setw(11) << analysis
                  << std::right;
    }

    file_path << 1 << '\n';

    // Describe the field and its scalar/vector/tensor components
    file_path << " -4  "
              << std::left
              << std::setw(8) << frd.name
              << std::right
              << std::setw(5) << static_cast<int>(frd.components.size())
              << std::setw(5) << 1
              << '\n';

    for (const FRDComponent& component : frd.components) {
        file_path << " -5  "
                  << std::left
                  << std::setw(8) << component.name
                  << std::right
                  << std::setw(5) << 1
                  << std::setw(5) << component.entity
                  << std::setw(5) << component.index_1
                  << std::setw(5) << component.index_2;

        if (component.derived) {
            file_path << std::setw(5)
                      << 1
                      << component.name;
        }

        file_path << '\n';
    }

    // Emit values in dense row order but address each row by its FRD node id
    constexpr Index values_per_line = 6;

    for (Index row = 0; row < field.rows; ++row) {
        for (Index component = 0;
             component < static_cast<Index>(frd.components.size());
             ++component) {
            if (component == 0) {
                file_path << " -1"
                          << std::setw(10)
                          << static_cast<long long>(
                                 node_ids[static_cast<std::size_t>(row)]);
            } else if (component % values_per_line == 0) {
                file_path << '\n'
                          << " -2"
                          << std::string(10, ' ');
            }

            write_float(file_path,
                field_value(field, row, component, frd.components[component]));
        }

        file_path << '\n';
    }

    file_path << " -3\n";
}

} // namespace writer
} // namespace io
} // namespace fem