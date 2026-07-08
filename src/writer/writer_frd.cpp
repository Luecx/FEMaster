#include "writer_frd.h"

#include "../model/element/element.h"
#include "../model/model_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <set>
#include <utility>

namespace fem {
namespace reader {

namespace {

std::string frd_name(const std::string& input,
                     std::size_t max_length,
                     const std::string& fallback) {
    std::string out;

    for (char c : input) {
        const unsigned char uc = static_cast<unsigned char>(c);

        if (std::isalnum(uc) || c == '_' || c == '-') {
            out.push_back(static_cast<char>(std::toupper(uc)));
        }

        if (out.size() == max_length) {
            break;
        }
    }

    return out.empty() ? fallback : out;
}

bool starts_with(const std::string& text,
                 const std::string& prefix) {
    return text.size() >= prefix.size()
        && text.compare(0, prefix.size(), prefix) == 0;
}

std::string dataset_name(const std::string& field_name) {
    const std::string name = frd_name(field_name, 64, "FIELD");

    if (name == "U" || name == "DISP" || name == "DISPLACEMENT") {
        return "DISP";
    }

    if (name == "RF" || name == "REACTION" || name == "REACTION_FORCE" ||
        name == "REACTION_FORCES") {
        return "RF";
    }

    if (name == "F" || name == "FORC" || name == "FORCE" || name == "FORCES" ||
        name == "LOAD" || name == "CLOAD" || name == "CLOADS" ||
        name == "EXTERNAL" || name == "EXTERNAL_FORCE" || name == "EXTERNAL_FORCES") {
        return "FORC";
    }

    if (name == "V" || name == "VELO" || name == "VELOCITY") {
        return "VELO";
    }

    if (name == "A" || name == "ACC" || name == "ACCELERATION") {
        return "ACC";
    }

    return frd_name(name, 8, "FIELD");
}

bool is_tensor_dataset(const std::string& name,
                       Index components) {
    return components == 6
        && (starts_with(name, "STRESS") || starts_with(name, "STRAIN"));
}

Index exported_components(const std::string& name,
                          Index field_components) {
    if (field_components <= 0) {
        return 0;
    }

    if (name == "DISP" && field_components >= 3) {
        return 4; // D1, D2, D3, ALL
    }

    if ((name == "RF" || name == "FORC" || name == "VELO" || name == "ACC") &&
        field_components >= 3) {
        return 3;
    }

    return std::min<Index>(field_components, static_cast<Index>(6));
}

std::string component_name(const std::string& name,
                           Index component,
                           Index components) {
    if (name == "DISP") {
        if (component == 0) return "D1";
        if (component == 1) return "D2";
        if (component == 2) return "D3";
        if (component == 3) return "ALL";
    }

    if (name == "RF") {
        if (component == 0) return "RF1";
        if (component == 1) return "RF2";
        if (component == 2) return "RF3";
    }

    if (name == "FORC") {
        if (component == 0) return "F1";
        if (component == 1) return "F2";
        if (component == 2) return "F3";
    }

    if (name == "VELO") {
        if (component == 0) return "V1";
        if (component == 1) return "V2";
        if (component == 2) return "V3";
    }

    if (name == "ACC") {
        if (component == 0) return "A1";
        if (component == 1) return "A2";
        if (component == 2) return "A3";
    }

    if (is_tensor_dataset(name, components)) {
        if (starts_with(name, "STRAIN")) {
            if (component == 0) return "EXX";
            if (component == 1) return "EYY";
            if (component == 2) return "EZZ";
            if (component == 3) return "EXY";
            if (component == 4) return "EYZ";
            if (component == 5) return "EZX";
        }

        if (component == 0) return "SXX";
        if (component == 1) return "SYY";
        if (component == 2) return "SZZ";
        if (component == 3) return "SXY";
        if (component == 4) return "SYZ";
        if (component == 5) return "SZX";
    }

    if (component == 0) return "C1";
    if (component == 1) return "C2";
    if (component == 2) return "C3";
    if (component == 3) return "C4";
    if (component == 4) return "C5";
    if (component == 5) return "C6";

    return "C";
}

void component_definition(const std::string& name,
                          Index component,
                          Index components,
                          int& entity,
                          int& index_1,
                          int& index_2,
                          bool& exists,
                          std::string& exists_name) {
    entity      = 1;
    index_1     = static_cast<int>(component + 1);
    index_2     = 0;
    exists      = false;
    exists_name = "";

    if (name == "DISP") {
        entity = 2;

        if (component == 3) {
            index_1     = 0;
            index_2     = 0;
            exists      = true;
            exists_name = "ALL";
        }

        return;
    }

    if (components == 3) {
        entity = 2;
        return;
    }

    if (is_tensor_dataset(name, components)) {
        entity = 4;

        if (component == 0) { index_1 = 1; index_2 = 1; return; }
        if (component == 1) { index_1 = 2; index_2 = 2; return; }
        if (component == 2) { index_1 = 3; index_2 = 3; return; }
        if (component == 3) { index_1 = 1; index_2 = 2; return; }
        if (component == 4) { index_1 = 2; index_2 = 3; return; }
        if (component == 5) { index_1 = 3; index_2 = 1; return; }
    }
}

double safe_value(double value) {
    return std::isfinite(value) ? value : 0.0;
}

double field_value(const model::Field& field,
                   Index row,
                   Index component,
                   const std::string& name) {
    if (name == "DISP" && component == 3) {
        const double u = safe_value(static_cast<double>(field(row, 0)));
        const double v = safe_value(static_cast<double>(field(row, 1)));
        const double w = safe_value(static_cast<double>(field(row, 2)));

        return std::sqrt(u * u + v * v + w * w);
    }

    return safe_value(static_cast<double>(field(row, component)));
}

void write_frd_float(std::ofstream& file_path,
                     double value) {
    file_path << std::right
              << std::uppercase
              << std::scientific
              << std::setprecision(5)
              << std::setw(12)
              << safe_value(value);
}

void write_result_value(std::ofstream& file_path,
                        double value) {
    file_path << std::right
              << std::uppercase
              << std::scientific
              << std::setprecision(5)
              << std::setw(11)
              << safe_value(value);
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

FrdWriter::FrdWriter(FrdWriter&& other) noexcept
    : file_path(std::move(other.file_path)),
      model_data_written(other.model_data_written),
      footer_written(other.footer_written),
      current_loadcase(other.current_loadcase),
      current_result_block(other.current_result_block),
      remap_zero_node_ids(other.remap_zero_node_ids),
      remap_zero_element_ids(other.remap_zero_element_ids),
      frd_node_ids(std::move(other.frd_node_ids)) {
    other.model_data_written   = false;
    other.footer_written       = true;
    other.current_result_block = 0;
}

FrdWriter& FrdWriter::operator=(FrdWriter&& other) noexcept {
    if (this != &other) {
        close();

        file_path              = std::move(other.file_path);
        model_data_written     = other.model_data_written;
        footer_written         = other.footer_written;
        current_loadcase       = other.current_loadcase;
        current_result_block   = other.current_result_block;
        remap_zero_node_ids    = other.remap_zero_node_ids;
        remap_zero_element_ids = other.remap_zero_element_ids;
        frd_node_ids           = std::move(other.frd_node_ids);

        other.model_data_written   = false;
        other.footer_written       = true;
        other.current_result_block = 0;
    }

    return *this;
}

void FrdWriter::open(const std::string& filename) {
    close();

    file_path.open(filename, std::ios::out | std::ios::trunc);

    logging::error(file_path.is_open(),
                   "FrdWriter: failed to open file: ", filename);

    model_data_written     = false;
    footer_written         = false;
    current_loadcase       = 1;
    current_result_block   = 0;
    remap_zero_node_ids    = false;
    remap_zero_element_ids = false;
    frd_node_ids.clear();

    write_header(filename);
}

void FrdWriter::close() {
    if (file_path.is_open()) {
        write_footer();
        file_path.close();
    }
}

bool FrdWriter::has_model_data() const {
    return model_data_written;
}

void FrdWriter::write_header(const std::string&) {
    file_path << "    1CFEMAST\n";
    file_path << "    1Ugenerated by FEMaster\n";
}

void FrdWriter::write_footer() {
    if (!footer_written) {
        file_path << " 9999\n";
        footer_written = true;
    }
}

void FrdWriter::add_loadcase(int id) {
    logging::error(file_path.is_open(),
                   "FrdWriter: cannot add loadcase: file is not open");

    current_loadcase = id > 0 ? id : 1;
}

void FrdWriter::write_model_data(const model::ModelData& model_data) {
    logging::error(file_path.is_open(),
                   "FrdWriter: cannot write model data: file is not open");

    if (model_data_written) {
        return;
    }

    logging::error(model_data.positions != nullptr,
                   "FrdWriter: model_data.positions is not initialized");

    frd_node_ids = collect_frd_node_ids(model_data);

    remap_zero_node_ids = std::find(frd_node_ids.begin(),
                                    frd_node_ids.end(),
                                    static_cast<ID>(0)) != frd_node_ids.end();

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
    logging::error(file_path.is_open(),
                   "FrdWriter: cannot write field '", field_name,
                   "': file is not open");

    if (!supports_field(field, field_name)) {
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

bool FrdWriter::supports_field(const model::Field& field,
                               const std::string&) const {
    return field.domain == model::FieldDomain::NODE;
}

std::vector<ID> FrdWriter::collect_frd_node_ids(const model::ModelData& model_data) const {
    std::set<ID> node_ids;

    for (const auto& element : model_data.elements) {
        if (!element || frd_element_type(*element) == 0) {
            continue;
        }

        const ID* nodes = element->nodes();

        for (Dim local = 0; local < element->n_nodes(); ++local) {
            node_ids.insert(nodes[static_cast<Index>(local)]);
        }
    }

    return std::vector<ID>(node_ids.begin(), node_ids.end());
}

std::size_t FrdWriter::count_supported_elements(const model::ModelData& model_data) const {
    std::size_t count = 0;

    for (const auto& element : model_data.elements) {
        if (element && frd_element_type(*element) != 0) {
            ++count;
        }
    }

    return count;
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

    logging::error(!frd_node_ids.empty(),
                   "FrdWriter: no supported nodes found for FRD output");

    file_path << "    2C"
              << std::string(18, ' ')
              << std::setw(12) << static_cast<long long>(frd_node_ids.size())
              << std::string(37, ' ')
              << 1
              << '\n';

    for (ID internal_node_id : frd_node_ids) {
        const Index row = static_cast<Index>(internal_node_id);

        logging::error(row >= 0 && row < positions.rows,
                       "FrdWriter: node id ", internal_node_id,
                       " is outside positions field rows");

        file_path << " -1"
                  << std::setw(10)
                  << static_cast<long long>(frd_node_number(internal_node_id));

        write_frd_float(file_path, static_cast<double>(positions(row, 0)));
        write_frd_float(file_path, static_cast<double>(positions(row, 1)));
        write_frd_float(file_path, static_cast<double>(positions(row, 2)));

        file_path << '\n';
    }

    file_path << " -3\n";
}

void FrdWriter::write_elements(const model::ModelData& model_data) {
    const std::size_t element_count = count_supported_elements(model_data);

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

        const ID* nodes = element->nodes();

        for (Dim local = 0; local < element->n_nodes(); local += static_cast<Dim>(10)) {
            file_path << " -2";

            const Dim n = std::min<Dim>(
                static_cast<Dim>(10),
                static_cast<Dim>(element->n_nodes() - local)
            );

            for (Dim k = 0; k < n; ++k) {
                const ID internal_node_id = nodes[static_cast<Index>(local + k)];

                file_path << std::setw(10)
                          << static_cast<long long>(frd_node_number(internal_node_id));
            }

            file_path << '\n';
        }
    }

    file_path << " -3\n";
}

void FrdWriter::write_nodal_field(const model::Field& field,
                                  const std::string& field_name) {
    logging::error(field.components > 0,
                   "FrdWriter: field '", field_name, "' has no components");

    logging::error(!frd_node_ids.empty(),
                   "FrdWriter: no FRD nodes available for field '", field_name, "'");

    const std::string name = dataset_name(field_name);
    const Index components = exported_components(name, field.components);

    logging::error(components > 0,
                   "FrdWriter: field '", field_name, "' has no exportable components");

    const int step  = 1;
    const int frame = current_loadcase > 0 ? current_loadcase : 1;
    const int block = ++current_result_block;

    file_path << "    1PSTEP"
              << std::setw(26) << block
              << std::setw(12) << frame
              << std::setw(12) << step
              << '\n';

    file_path << "  100CL  "
              << std::setw(3) << 100 + block
              << ' ';

    write_result_value(file_path, static_cast<double>(frame));

    file_path << std::setw(12) << static_cast<long long>(frd_node_ids.size())
              << std::string(21, ' ')
              << 0
              << std::setw(5) << block
              << std::string(11, ' ')
              << 1
              << '\n';

    file_path << " -4  "
              << std::left << std::setw(8) << name
              << std::right
              << std::setw(5) << static_cast<int>(components)
              << std::setw(5) << 1
              << '\n';

    for (Index component = 0; component < components; ++component) {
        int entity;
        int index_1;
        int index_2;
        bool exists;
        std::string exists_name;

        component_definition(name,
                             component,
                             components,
                             entity,
                             index_1,
                             index_2,
                             exists,
                             exists_name);

        const std::string label = component_name(name, component, components);

        file_path << " -5  "
                  << std::left << std::setw(8) << label
                  << std::right
                  << std::setw(5) << 1
                  << std::setw(5) << entity
                  << std::setw(5) << index_1
                  << std::setw(5) << index_2;

        if (exists) {
            file_path << std::setw(5) << 1 << exists_name;
        }

        file_path << '\n';
    }

    for (ID internal_node_id : frd_node_ids) {
        const Index row = static_cast<Index>(internal_node_id);

        logging::error(row >= 0 && row < field.rows,
                       "FrdWriter: node id ", internal_node_id,
                       " is outside rows of field '", field_name, "'");

        file_path << " -1"
                  << std::setw(10)
                  << static_cast<long long>(frd_node_number(internal_node_id));

        for (Index component = 0; component < components; ++component) {
            write_frd_float(file_path, field_value(field, row, component, name));
        }

        file_path << '\n';
    }

    file_path << " -3\n";
}

int FrdWriter::frd_element_type(const model::ElementInterface& element) const {
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

} // namespace reader
} // namespace fem