/**
 * @file writer_frd.cpp
 * @brief Implements the ASCII CalculiX/CGX FRD result writer.
 *
 * The writer emits compiled FEMaster nodes, supported finite elements and
 * NODE-domain structural or thermal result fields in CalculiX FRD syntax.
 *
 * FEMaster solver topology uses dense node identifiers. When model data is
 * written, the semantic reverse mapping retained by `ModelData` is evaluated
 * once to construct the external FRD identifier
 *
 *     n_FRD = 10^8 i_instance + n_local.
 *
 * The resulting dense-to-FRD mapping is cached locally in `FrdWriter`. Node
 * coordinates, element connectivity and every subsequent nodal result block
 * therefore use identical external identifiers without repeated Instance
 * lookups.
 *
 * Large mesh and result blocks are formatted into independent thread-local
 * string buffers. The completed buffers are written sequentially in their
 * original order so parallel formatting does not alter the FRD record order.
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

#include "../../core/config.h"
#include "../../model/element/element.h"
#include "../../model/element/element_structural.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <string>
#include <string_view>
#include <system_error>
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

/**
 * @brief Defines one supported FRD result field.
 *
 * Aliases contain the normalized FEMaster names accepted for this field.
 * Components describe the scalar, vector and tensor metadata emitted through
 * FRD -5 records.
 */
struct FRDField {
    std::vector<std::string> aliases;
    std::string name;
    std::vector<FRDComponent> components;
};


/**
 * Appends a left-aligned string using FRD fixed-width formatting.
 */
void append_string_left(std::string& buffer,
                        std::string_view value,
                        std::size_t width) {
    buffer.append(value);

    if (value.size() < width) {
        buffer.append(width - value.size(), ' ');
    }
}

/**
 * Appends a right-aligned integer using locale-independent conversion.
 */
template<typename T>
void append_integer(std::string& buffer,
                    T value,
                    std::size_t width = 0) {
    char tmp[32];

    const auto [end, ec] = std::to_chars(tmp, tmp + sizeof(tmp), value);

    logging::error(ec == std::errc{},
        "FrdWriter: failed to format integer");

    const std::size_t length = static_cast<std::size_t>(end - tmp);

    if (length < width) {
        buffer.append(width - length, ' ');
    }

    buffer.append(tmp, length);
}

/**
 * Appends one FRD floating-point value in uppercase scientific notation.
 *
 * FRD numerical records use five digits after the decimal point. Non-finite
 * values are replaced by zero before formatting, preserving the previous writer
 * behavior.
 */
void append_float(std::string& buffer,
                  Precision value,
                  std::size_t width = 12) {
    char tmp[64];

    if (!std::isfinite(value)) {
        value = Precision(0);
    }

    const auto [end, ec] = std::to_chars(
        tmp,
        tmp + sizeof(tmp),
        value,
        std::chars_format::scientific,
        5
    );

    logging::error(ec == std::errc{},
        "FrdWriter: failed to format floating-point value");

    const std::size_t length = static_cast<std::size_t>(end - tmp);

    if (length < width) {
        buffer.append(width - length, ' ');
    }

    for (char* ptr = tmp; ptr != end; ++ptr) {
        buffer.push_back(*ptr == 'e' ? 'E' : *ptr);
    }
}

/**
 * Selects the number of formatting workers for one output block.
 *
 * Approximately one worker is enabled per 4096 rows or elements. The count is
 * bounded by the configured FEMaster thread count and by 32 writer threads.
 */
ID output_thread_count(Index work_items) {
#ifndef _OPENMP
    (void) work_items;
    return 1;
#else
    const Index blocks = std::max<Index>(1, (work_items + 4095) / 4096);
    const ID useful_threads = static_cast<ID>(std::min<Index>(blocks, 32));

    return std::max<ID>(1, std::min<ID>(global_config.max_threads, useful_threads));
#endif
}

/**
 * Creates one independent text buffer for every formatting worker.
 */
std::vector<std::string> create_thread_buffers(ID num_threads) {
    std::vector<std::string> buffers(static_cast<std::size_t>(num_threads));

    for (auto& buffer : buffers) {
        buffer.reserve(4 * 1024 * 1024);
    }

    return buffers;
}

/**
 * Writes completed worker buffers in their original block order.
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
            {"TEMPERATURE", "NDTEMP"}, "NDTEMP",
            {
                FRDComponent::scalar("T", 1)
            }
        },
        {
            {"HEATFLUX", "HFL", "FLUX"}, "FLUX",
            {
                FRDComponent::vector("HFL1", 1),
                FRDComponent::vector("HFL2", 2),
                FRDComponent::vector("HFL3", 3),
                FRDComponent::vcnorm("ALL")
            }
        },
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
            {"DISPLACEMENTREAL"}, "DISPREAL",
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
            {"DISPLACEMENTIMAG"}, "DISPIMAG",
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
            {"STRESSREAL"}, "STRREAL",
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
            {"STRESSIMAG"}, "STRIMAG",
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
            {"STRAINREAL"}, "STRNREAL",
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
            {"STRAINIMAG"}, "STRNIMAG",
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
 * Appends the FRD node-block header.
 */
void append_node_block_header(std::string& buffer, Index node_count) {
    buffer.append("    2C");
    buffer.append(18, ' ');
    append_integer(buffer, node_count, 12);
    buffer.append(37, ' ');
    buffer.push_back('1');
    buffer.push_back('\n');
}

/**
 * Appends one -1 node-coordinate record.
 */
void append_node_line(std::string& buffer,
                      ID node_id,
                      Precision x,
                      Precision y,
                      Precision z) {
    buffer.append(" -1");
    append_integer(buffer, node_id, 10);
    append_float(buffer, x);
    append_float(buffer, y);
    append_float(buffer, z);
    buffer.push_back('\n');
}

/**
 * Appends the FRD element-block header.
 */
void append_element_block_header(std::string& buffer,
                                 std::size_t element_count) {
    buffer.append("    3C");
    buffer.append(18, ' ');
    append_integer(buffer, element_count, 12);
    buffer.append(37, ' ');
    buffer.push_back('1');
    buffer.push_back('\n');
}

/**
 * Appends one -1 element-definition record.
 */
void append_element_header_line(std::string& buffer,
                                ID element_id,
                                int type) {
    buffer.append(" -1");
    append_integer(buffer, element_id, 10);
    append_integer(buffer, type, 5);
    append_integer(buffer, 0, 5);
    append_integer(buffer, 0, 5);
    buffer.push_back('\n');
}

/**
 * Appends the -2 connectivity records belonging to one element.
 *
 * FRD stores at most ten node identifiers on one connectivity line. Additional
 * nodes continue on a new -2 record.
 */
void append_element_connectivity(std::string& buffer,
                                 const std::vector<ID>& connectivity,
                                 const std::vector<ID>& node_ids,
                                 ID element_id) {
    for (std::size_t local = 0; local < connectivity.size(); ++local) {
        const ID node_id = connectivity[local];

        logging::error(static_cast<std::size_t>(node_id) < node_ids.size(),
            "FrdWriter: element ", element_id,
            " references node ", node_id,
            " outside the FRD node mapping");

        if (local % 10 == 0) {
            if (local != 0) {
                buffer.push_back('\n');
            }

            buffer.append(" -2");
        }

        append_integer(
            buffer,
            node_ids[static_cast<std::size_t>(node_id)],
            10
        );
    }

    buffer.push_back('\n');
}

/**
 * Appends one 1PSTEP result-step metadata record.
 */
void append_result_step_line(std::string& buffer,
                             int block,
                             int frame,
                             int step) {
    buffer.append("    1PSTEP");
    append_integer(buffer, block, 26);
    append_integer(buffer, frame, 12);
    append_integer(buffer, step, 12);
    buffer.push_back('\n');
}

/**
 * Appends one 1PMODE metadata record.
 */
void append_result_mode_line(std::string& buffer, int frame) {
    buffer.append("    1PMODE");
    append_integer(buffer, frame, 26);
    buffer.push_back('\n');
}

/**
 * Appends the 100CL result-control record.
 */
void append_result_control_line(std::string& buffer,
                                int global_frame,
                                Precision value,
                                Index rows,
                                int ictype,
                                const char* analysis) {
    buffer.append("  100CL  ");
    append_integer(buffer, 100 + global_frame, 3);
    buffer.push_back(' ');

    append_float(buffer, value, 11);
    append_integer(buffer, rows, 12);

    buffer.append(21, ' ');
    append_integer(buffer, ictype);
    append_integer(buffer, global_frame, 5);

    if (analysis[0] == '\0') {
        buffer.append(11, ' ');
    } else {
        append_string_left(buffer, analysis, 11);
    }

    buffer.push_back('1');
    buffer.push_back('\n');
}

/**
 * Appends the -4 field-definition record.
 */
void append_field_definition_line(std::string& buffer,
                                  const FRDField& field) {
    buffer.append(" -4  ");
    append_string_left(buffer, field.name, 8);
    append_integer(buffer, static_cast<int>(field.components.size()), 5);
    append_integer(buffer, 1, 5);
    buffer.push_back('\n');
}

/**
 * Appends one -5 component-definition record.
 */
void append_component_definition_line(std::string& buffer,
                                      const FRDComponent& component) {
    buffer.append(" -5  ");
    append_string_left(buffer, component.name, 8);
    append_integer(buffer, 1, 5);
    append_integer(buffer, component.entity, 5);
    append_integer(buffer, component.index_1, 5);
    append_integer(buffer, component.index_2, 5);

    if (component.derived) {
        append_integer(buffer, 1, 5);
        buffer.append(component.name);
    }

    buffer.push_back('\n');
}

/**
 * Appends all -1/-2 numerical result records for one node.
 *
 * The first result line starts with the external node identifier. Additional
 * components continue on -2 records after every six values.
 */
void append_result_node_lines(std::string& buffer,
                              ID node_id,
                              const model::Field& field,
                              Index row,
                              const FRDField& frd) {
    constexpr Index values_per_line = 6;

    for (Index component = 0;
         component < static_cast<Index>(frd.components.size());
         ++component) {

        if (component == 0) {
            buffer.append(" -1");
            append_integer(buffer, node_id, 10);
        } else if (component % values_per_line == 0) {
            buffer.push_back('\n');
            buffer.append(" -2");
            buffer.append(10, ' ');
        }

        append_float(
            buffer,
            field_value(field, row, component, frd.components[component])
        );
    }

    buffer.push_back('\n');
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
            "FrdWriter: local node id ", local_id,
            " exceeds the supported instance range");

        node_ids[i] =
            static_cast<ID>(100000000) * static_cast<ID>(instance->instance_id)
            + local_id;
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
 * Coordinate rows are formatted in contiguous blocks by independent workers.
 * The buffers are subsequently written in block order to preserve dense node
 * ordering.
 *
 * @param model_data Compiled model providing the nodal position field.
 */
void FrdWriter::write_nodes(const model::ModelData& model_data) {
    const auto& positions = *model_data.positions;

    logging::error(node_ids.size() == static_cast<std::size_t>(positions.rows),
        "FrdWriter: node id mapping does not match positions field");

    // Write the block header before the parallel node data
    std::string header;
    header.reserve(128);
    append_node_block_header(header, positions.rows);
    file_path.write(header.data(), static_cast<std::streamsize>(header.size()));

    // Format contiguous node ranges independently
    const ID num_threads = output_thread_count(positions.rows);
    auto thread_buffers = create_thread_buffers(num_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
#endif
    for (ID thread = 0; thread < num_threads; ++thread) {
        auto& buffer = thread_buffers[static_cast<std::size_t>(thread)];

        const Index begin = positions.rows * thread       / num_threads;
        const Index end   = positions.rows * (thread + 1) / num_threads;

        for (Index row = begin; row < end; ++row) {
            append_node_line(
                buffer,
                node_ids[static_cast<std::size_t>(row)],
                positions(row, 0),
                positions(row, 1),
                positions(row, 2)
            );
        }
    }

    // Commit node ranges in their original dense order
    write_thread_buffers(file_path, thread_buffers);

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
 * Node-only models omit the connectivity block while retaining their nodal
 * geometry and result fields.
 *
 * Each worker formats complete elements so FRD connectivity records never cross
 * thread-buffer boundaries.
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

    // Keep node-only feature models valid without emitting an empty element block
    if (element_count == 0) {
        return;
    }

    // Write the element block header before the parallel element data
    std::string header;
    header.reserve(128);
    append_element_block_header(header, element_count);
    file_path.write(header.data(), static_cast<std::streamsize>(header.size()));

    // Divide the compiled element range into contiguous worker blocks
    const Index compiled_count = static_cast<Index>(model_data.elements.size());
    const ID num_threads = output_thread_count(compiled_count);
    auto thread_buffers = create_thread_buffers(num_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
#endif
    for (ID thread = 0; thread < num_threads; ++thread) {
        auto& buffer = thread_buffers[static_cast<std::size_t>(thread)];

        const Index begin = compiled_count * thread       / num_threads;
        const Index end   = compiled_count * (thread + 1) / num_threads;

        for (Index elem_row = begin; elem_row < end; ++elem_row) {
            const auto& element =
                model_data.elements[static_cast<std::size_t>(elem_row)];

            if (!element) {
                continue;
            }

            const int type = frd_element_type(*element);

            if (type == 0) {
                continue;
            }

            append_element_header_line(buffer, element->elem_id, type);

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

            append_element_connectivity(
                buffer,
                connectivity,
                node_ids,
                element->elem_id
            );
        }
    }

    // Commit complete element ranges in compiled order
    write_thread_buffers(file_path, thread_buffers);

    file_path << " -3\n";
}

/**
 * Writes one NODE-domain result block using the cached FRD node numbering.
 *
 * Field rows correspond directly to dense global node identifiers. Each row is
 * therefore translated by indexing the same `node_ids` vector used for mesh
 * connectivity, guaranteeing consistent numbering across geometry and results.
 *
 * Result metadata is small and remains serial. The potentially large -1/-2
 * numerical record block is formatted in parallel using contiguous node ranges.
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

    // Build all compact result metadata records in one local buffer
    std::string metadata;
    metadata.reserve(1024);

    append_result_step_line(metadata, block, frame, step);

    if (current_step_type == WriterStepType::Eigenfrequency) {
        append_result_mode_line(metadata, frame);
    }

    append_result_control_line(
        metadata,
        global_frame,
        value,
        field.rows,
        ictype,
        analysis
    );

    append_field_definition_line(metadata, frd);

    for (const FRDComponent& component : frd.components) {
        append_component_definition_line(metadata, component);
    }

    file_path.write(
        metadata.data(),
        static_cast<std::streamsize>(metadata.size())
    );

    // Format contiguous node-result ranges independently
    const ID num_threads = output_thread_count(field.rows);
    auto thread_buffers = create_thread_buffers(num_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
#endif
    for (ID thread = 0; thread < num_threads; ++thread) {
        auto& buffer = thread_buffers[static_cast<std::size_t>(thread)];

        const Index begin = field.rows * thread       / num_threads;
        const Index end   = field.rows * (thread + 1) / num_threads;

        for (Index row = begin; row < end; ++row) {
            append_result_node_lines(
                buffer,
                node_ids[static_cast<std::size_t>(row)],
                field,
                row,
                frd
            );
        }
    }

    // Preserve dense node order in the completed result block
    write_thread_buffers(file_path, thread_buffers);

    file_path << " -3\n";
}

} // namespace writer
} // namespace io
} // namespace fem
