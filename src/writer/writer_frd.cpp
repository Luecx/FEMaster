#include "writer_frd.h"

#include "../model/element/element.h"
#include "../model/model_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <set>

namespace fem {
namespace reader {

namespace {

struct FRDComponent {
    std::string name;
    int entity      = 1; // scalar = 1, vector = 2, tensor = 4
    int index_1     = 0; // first index, e.g. for vectors or tensors
    int index_2     = 0; // only relevant for tensors
    std::string exists_name;

    static FRDComponent scalar(const std::string& name, int index   ) { return {name, 1, index, 0, ""};}
    static FRDComponent vector(const std::string& name, int index   ) { return {name, 2, index, 0, ""};}
    static FRDComponent tensor(const std::string& name, int i, int j) { return {name, 4, i, j, ""};}
    static FRDComponent vcnorm(const std::string& name)               { return {name, 2, 0, 0, name};}
};

int frd_frame(const std::string& s) {
    std::string out;
    for (unsigned char c : s) {
        if (std::isdigit(c)) out += static_cast<char>(c);
    }
    return out.empty() ? 1 : std::stoi(out);
}


int frd_element_type(const model::ElementInterface& element) {
    const Dim dim     = element.dimensions();
    const Dim n_nodes = element.n_nodes();

    if (dim == 1) {
        if (n_nodes == 2) return 11;
        if (n_nodes == 3) return 12;
    }

    if (dim == 2) {
        if (n_nodes == 3) return 7;
        if (n_nodes == 6) return 8;
        if (n_nodes == 4) return 9;
        if (n_nodes == 8) return 10;
    }

    if (dim == 3) {
        if (n_nodes == 2)  return 11;
        if (n_nodes == 3)  return 12;
        if (n_nodes == 4)  return 3;
        if (n_nodes == 10) return 6;
        if (n_nodes == 8)  return 1;
        if (n_nodes == 20) return 4;
        if (n_nodes == 6)  return 2;
        if (n_nodes == 15) return 5;
    }

    return 0;
}

std::string frd_name(const std::string& s) {
    std::string name;
    for (unsigned char c : s) {
        if (std::isalpha(c)) {
            name += static_cast<char>(std::toupper(c));
        }
    }

    if (name == "DISPLACEMENT"      ) return "DISP";
    if (name == "MODESHAPE"         ) return "DISP";
    if (name == "BUCKLINGMODE"      ) return "DISP";

    if (name == "VELOCITY"          ) return "V";
    if (name == "ACCELERATION"      ) return "A";

    if (name == "STRESS"            ) return "STRESS";
    if (name == "STRAIN"            ) return "STRAIN";
    if (name == "STRESSTOP"         ) return "S_TOP";
    if (name == "STRESSBOT"         ) return "S_BOT";

    if (name == "REACTIONFORCES"    ) return "RF";
    if (name == "EXTERNALFORCES"    ) return "CF";
    if (name == "INTERNALFORCES"    ) return "NFORC";

    if (name == "LOCALSECTIONFORCES") return "SF";
    if (name == "SHELLRESULTANTS"   ) return "SHR";
    if (name == "SHEARFLOW"         ) return "SHEAR";

    logging::error(false, "FrdWriter: unsupported NODE field name: ", s);
    return "FIELD";
}

std::vector<FRDComponent> vector_components(const std::string& prefix,
                                            Index field_components,
                                            bool add_norm = false) {
    std::vector<FRDComponent> components;

    for (Index i = 0; i < std::min<Index>(field_components, 3); ++i) {
        components.push_back(FRDComponent::vector(prefix + std::to_string(i + 1), static_cast<int>(i + 1)));
    }

    if (add_norm && field_components >= 3) {
        components.push_back(FRDComponent::vcnorm("ALL"));
    }

    return components;
}

std::vector<FRDComponent> tensor_components(const std::string& prefix,
                                            Index field_components) {
    logging::error(field_components >= 6,
                   "FrdWriter: tensor field requires at least 6 components");

    return {
        FRDComponent::tensor(prefix + "XX", 1, 1),
        FRDComponent::tensor(prefix + "YY", 2, 2),
        FRDComponent::tensor(prefix + "ZZ", 3, 3),
        FRDComponent::tensor(prefix + "YZ", 2, 3),
        FRDComponent::tensor(prefix + "ZX", 3, 1),
        FRDComponent::tensor(prefix + "XY", 1, 2)
    };
}

std::vector<FRDComponent> scalar_components(const std::string& prefix,
                                            Index field_components) {
    std::vector<FRDComponent> components;

    for (Index i = 0; i < field_components; ++i) {
        components.push_back(FRDComponent::scalar(prefix + std::to_string(i + 1), static_cast<int>(i + 1)));
    }

    return components;
}

std::vector<FRDComponent> components_for(const std::string& name,
                                         Index field_components) {
    logging::error(field_components > 0,
                   "FrdWriter: field has no components");

    if (name == "DISP"  ) return vector_components("D" , field_components, true);
    if (name == "V"     ) return vector_components("V" , field_components);
    if (name == "A"     ) return vector_components("A" , field_components);
    if (name == "RF"    ) return vector_components("RF", field_components);
    if (name == "CF"    ) return vector_components("CF", field_components);
    if (name == "NFORC" ) return vector_components("NF", field_components);

    if (name == "STRESS") return tensor_components("S", field_components);
    if (name == "STRAIN") return tensor_components("E", field_components);
    if (name == "S_TOP" ) return tensor_components("S", field_components);
    if (name == "S_BOT" ) return tensor_components("S", field_components);

    if (name == "SF"    ) return scalar_components("SF"   , field_components);
    if (name == "SHR"   ) return scalar_components("SHR"  , field_components);
    if (name == "SHEAR" ) return scalar_components("SHEAR", field_components);

    logging::error(false, "FrdWriter: unsupported FRD result name: ", name);
    return {};
}

Precision field_value(const model::Field& field,
                      Index row,
                      Index component,
                      const std::string& name) {
    if (name == "DISP" && component == 3) {
        const Precision u = std::isfinite(field(row, 0)) ? field(row, 0) : 0.0;
        const Precision v = std::isfinite(field(row, 1)) ? field(row, 1) : 0.0;
        const Precision w = std::isfinite(field(row, 2)) ? field(row, 2) : 0.0;

        return std::sqrt(u * u + v * v + w * w);
    }

    return std::isfinite(field(row, component)) ? field(row, component) : 0;
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

void FrdWriter::add_loadcase(int id) {
    logging::error(file_path.is_open(), "FrdWriter: cannot add loadcase: file is not open");

    current_step = id > 0 ? id : 1;
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
                            const model::ModelData* model_data) {
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

    write_nodal_field(field, field_name);
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

    logging::error(!frd_node_ids.empty(), "FrdWriter: no supported nodes found for FRD output");

    file_path << "    2C"
              << std::string(18, ' ')
              << std::setw(12) << static_cast<long long>(frd_node_ids.size())
              << std::string(37, ' ')
              << 1
              << '\n';

    for (ID node_id : frd_node_ids) {
        const Index row = static_cast<Index>(node_id);

        logging::error(row >= 0 && row < positions.rows,
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
              << std::setw(12) << static_cast<long long>(element_count)
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

        Dim local = 0;

        for (ID node_id : *element) {
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
                                  const std::string& field_name) {
    logging::error(field.components > 0 , "FrdWriter: field '", field_name, "' has no components");
    logging::error(!frd_node_ids.empty(), "FrdWriter: no FRD nodes available for field '", field_name, "'");

    const std::string name = frd_name(field_name);
    const int frame        = frd_frame(field_name);
    const int step         = current_step > 0 ? current_step : 1;
    const int block        = ++current_result_block;
    const auto components  = components_for(name, field.components);

    file_path << "    1PSTEP" << std::setw(26) <<       block << std::setw(12) << frame << std::setw(12) << step << '\n';
    file_path << "  100CL  "  << std::setw(3 ) << 100 + block << ' ';

    write_float(file_path, static_cast<Precision>(block), 11);

    file_path << std::setw  (12)      << frd_node_ids.size()
              << std::string(21, ' ') << 0
              << std::setw  (5)       << block
              << std::string(11, ' ') << 1
              << '\n';

    file_path << " -4  "
              << std::left << std::setw(8) << name
              << std::right
              << std::setw(5) << static_cast<int>(components.size())
              << std::setw(5) << 1
              << '\n';

    for (const FRDComponent& component : components) {
        file_path << " -5  "
                  << std::left << std::setw(8) << component.name
                  << std::right
                  << std::setw(5) << 1
                  << std::setw(5) << component.entity
                  << std::setw(5) << component.index_1
                  << std::setw(5) << component.index_2;

        if (!component.exists_name.empty()) {
            file_path << std::setw(5) << 1 << component.exists_name;
        }

        file_path << '\n';
    }

    for (ID node_id : frd_node_ids) {
        const Index row = static_cast<Index>(node_id);

        logging::error(row >= 0 && row < field.rows,
                       "FrdWriter: node id ", node_id,
                       " is outside rows of field '", field_name, "'");

        file_path << " -1" << std::setw(10) << frd_node_number(node_id);

        for (Index component = 0; component < static_cast<Index>(components.size()); ++component) {
            write_float(file_path, field_value(field, row, component, name));
        }

        file_path << '\n';
    }

    file_path << " -3\n";
}


} // namespace reader
} // namespace fem