#include "writer_frd.h"

#include "../model/element/element.h"
#include "../model/element/element_structural.h"
#include "../model/model_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <set>
#include <vector>

namespace fem {
namespace reader {

namespace {

// a single component, e.g. the stress.Sxy component.
struct FRDComponent {
    std::string name;
    int entity   = 1; // scalar = 1, vector = 2, tensor = 4
    int index_1  = 0; // first index, e.g. for vectors or tensors
    int index_2  = 0; // only relevant for tensors
    bool derived = false;

    static FRDComponent scalar(const std::string& name, int index   ) { return {name, 1, index, 0, false};}
    static FRDComponent vector(const std::string& name, int index   ) { return {name, 2, index, 0, false};}
    static FRDComponent tensor(const std::string& name, int i, int j) { return {name, 4, i, j, false};}
    static FRDComponent vcnorm(const std::string& name)               { return {name, 2, 0, 0, true };}
};

struct FRDField {
    std::vector<std::string> aliases;
    std::string name;
    std::vector<FRDComponent> components;
};

int frd_frame(const std::string& s) {
    std::string out;

    for (unsigned char c : s) {
        if (std::isdigit(c)) {
            out += static_cast<char>(c);
        }
    }

    return out.empty() ? 1 : std::stoi(out);
}

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

int frd_increment_type(WriterStepType step_type) {
    switch (step_type) {
        case WriterStepType::Dynamic:        return 1;
        case WriterStepType::Eigenfrequency: return 2;
        case WriterStepType::Buckling:       return 4;
        case WriterStepType::Static:         return 0;
    }

    return 0;
}

const char* frd_analysis_name(WriterStepType step_type) {
    switch (step_type) {
        case WriterStepType::Dynamic:        return "DYNAMIC";
        case WriterStepType::Eigenfrequency: return "MODAL";
        case WriterStepType::Buckling:       return "BUCKLING";
        case WriterStepType::Static:         return "";
    }

    return "";
}

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
                   "FrdWriter: unsupported NODE field name: ",
                   field_name);

    return definitions.front();
}

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

FrdWriter::FrdWriter(const std::string& filename) {
    if (!filename.empty()) {
        open(filename);
    }
}

FrdWriter::~FrdWriter() {
    close();
}

void FrdWriter::open(const std::string& filename) {
    close();

    file_path.open(filename, std::ios::out | std::ios::trunc);

    logging::error(file_path.is_open(),
                   "FrdWriter: failed to open file: ", filename);

    model_data_written     = false;
    current_step           = 1;
    current_result_block   = 0;
    current_global_frame   = 0;
    last_frame_step        = -1;
    last_frame_id          = -1;
    current_step_type      = WriterStepType::Static;
    remap_zero_node_ids    = false;
    remap_zero_element_ids = false;
    frd_node_ids.clear();

    file_path << "    1CFEMAST\n";
    file_path << "    1Ugenerated by FEMaster\n";
}

void FrdWriter::close() {
    if (file_path.is_open()) {
        file_path << " 9999\n";
        file_path.close();
    }
}

void FrdWriter::add_loadcase(int id, WriterStepType step_type) {
    logging::error(file_path.is_open(),
                   "FrdWriter: cannot add loadcase: file is not open");

    current_step      = id > 0 ? id : 1;
    current_step_type = step_type;
    last_frame_step   = -1;
    last_frame_id     = -1;
}

void FrdWriter::write_model_data(const model::ModelData& model_data) {
    logging::error(file_path.is_open()            , "FrdWriter: cannot write model data: file is not open");
    logging::error(model_data.positions != nullptr, "FrdWriter: model_data.positions is not initialized");

    if (model_data_written) {
        return;
    }

    remap_zero_element_ids = false;

    for (const auto& element : model_data.elements) {
        if (element &&
            frd_element_type(*element) != 0 &&
            element->elem_id == static_cast<ID>(0)) {
            remap_zero_element_ids = true;
            break;
        }
    }

    write_nodes(model_data);
    write_elements(model_data);

    model_data_written = true;
}

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

ID FrdWriter::frd_node_number(ID internal_node_id) const {
    return remap_zero_node_ids
        ? internal_node_id + static_cast<ID>(1)
        : internal_node_id;
}

ID FrdWriter::frd_element_number(ID internal_element_id) const {
    return remap_zero_element_ids
        ? internal_element_id + static_cast<ID>(1)
        : internal_element_id;
}

void FrdWriter::write_nodes(const model::ModelData& model_data) {
    const auto& positions = *model_data.positions;

    std::set<ID> node_ids;

    for (const auto& element : model_data.elements) {
        if (!element || frd_element_type(*element) == 0) {
            continue;
        }

        for (ID node_id : *element) {
            node_ids.insert(node_id);
        }
    }

    frd_node_ids.assign(node_ids.begin(), node_ids.end());

    remap_zero_node_ids = std::find(frd_node_ids.begin(),
                                    frd_node_ids.end(),
                                    static_cast<ID>(0)) != frd_node_ids.end();

    logging::error(!frd_node_ids.empty(),
                   "FrdWriter: no supported nodes found for FRD output");

    file_path << "    2C"
              << std::string(18, ' ')
              << std::setw(12)
              << static_cast<long long>(frd_node_ids.size())
              << std::string(37, ' ')
              << 1
              << '\n';

    for (ID node_id : frd_node_ids) {
        const Index row = static_cast<Index>(node_id);

        logging::error(row < positions.rows,
                       "FrdWriter: node id ", node_id,
                       " is outside positions field rows");

        file_path << " -1"
                  << std::setw(10)
                  << static_cast<long long>(frd_node_number(node_id));

        write_float(file_path, positions(row, 0));
        write_float(file_path, positions(row, 1));
        write_float(file_path, positions(row, 2));

        file_path << '\n';
    }

    file_path << " -3\n";
}

void FrdWriter::write_elements(const model::ModelData& model_data) {
    std::size_t element_count = 0;

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
                  << static_cast<long long>(frd_element_number(element->elem_id))
                  << std::setw(5) << type
                  << std::setw(5) << 0
                  << std::setw(5) << 0
                  << '\n';

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

        Dim local = 0;

        for (ID node_id : connectivity) {
            if (local % static_cast<Dim>(10) == 0) {
                if (local != 0) {
                    file_path << '\n';
                }

                file_path << " -2";
            }

            file_path << std::setw(10)
                      << static_cast<long long>(frd_node_number(node_id));

            ++local;
        }

        file_path << '\n';
    }

    file_path << " -3\n";
}
void FrdWriter::write_nodal_field(const model::Field& field,
                                  const std::string& field_name,
                                  Precision frame_value) {
    logging::error(!frd_node_ids.empty(),
                   "FrdWriter: no FRD nodes available for field '",
                   field_name, "'");

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

    file_path << std::setw(12) << frd_node_ids.size()
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

    constexpr Index values_per_line = 6;

    for (ID node_id : frd_node_ids) {
        const Index row = static_cast<Index>(node_id);

        logging::error(row < field.rows,
                       "FrdWriter: node id ", node_id,
                       " is outside rows of field '",
                       field_name, "'");

        for (Index component = 0;
             component < static_cast<Index>(frd.components.size());
             ++component) {
            if (component == 0) {
                file_path << " -1"
                          << std::setw(10)
                          << frd_node_number(node_id);
            } else if (component % values_per_line == 0) {
                file_path << '\n'
                          << " -2"
                          << std::string(10, ' ');
            }

            write_float(file_path, field_value(field, row, component, frd.components[component]));
             }

        file_path << '\n';
    }

    file_path << " -3\n";
}

} // namespace reader
} // namespace fem
