/**
 * @file register_imperfection.cpp
 * @brief Applies nodal geometric imperfections to compiled reference geometry.
 *
 * Direct nodal offsets and scaled displacement frames from FEMaster RES or FRD
 * results are accumulated in the global Cartesian configuration. The reference
 * and initial current coordinate fields are updated together, before reference
 * shell normals and other geometry-dependent analysis data are initialized.
 *
 * RES frames are identified by distinct displacement-field suffixes within
 * each LC block. FRD frames use explicit 1PSTEP metadata. Only the first three displacement
 * components contribute; rotations are not changes in nodal position.
 *
 * @author Finn Eggers
 * @date 04.11.2026
 */

#include "register_functions.h"
#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../../model/model.h"
#include "../../../core/logging.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem::io::reader::commands {

namespace {

/**
 * Reads requested one-based displacement frames from a FEMaster RES file.
 *
 * RES stores an LC number and individual named fields, but no explicit frame
 * record. Suffixed fields such as BUCKLING_MODE_1 and DISPLACEMENT_1 identify
 * the same written frame when their suffix matches. The first distinct
 * displacement-frame suffix encountered is frame 1, the next is frame 2, etc.
 * Unsuffixed displacement fields each define a separate frame.
 *
 * Only NODE-domain fields named DISPLACEMENT, BUCKLING_MODE or MODE_SHAPE
 * (with or without a numerical suffix) are considered. If multiple aliases
 * are present for one frame, the first is used once.
 *
 * @param filename Source RES path.
 * @param step LC identifier to read.
 * @param scales Requested one-based frame indices and scaling factors.
 * @param model Current compiled model resolving semantic node identifiers.
 * @param offsets Accumulated global Cartesian nodal imperfections.
 */
void read_res(const std::string&           filename,
              int                          step,
              const std::map<int, Precision>& scales,
              model::Model&                model,
              std::vector<Vec3>&           offsets) {
    // Open and validate the requested result source
    std::ifstream file(filename);

    logging::error(file.is_open(),
        "IMPERFECTION: cannot open RES file ", filename);

    // Track written frame identities and requested displacement blocks
    std::string               line;
    std::map<std::string, int> numbered_frames;
    std::set<int>             processed_frames;
    std::set<int>             found_frames;

    int  current_step   = -1;
    int  frame_count    = 0;
    int  selected_frame = 0;
    bool read_field     = false;

    while (std::getline(file, line)) {
        // A new LC resets frame numbering and field selection
        if (line.rfind("LC ", 0) == 0) {
            std::istringstream stream(line.substr(3));
            stream >> current_step;

            logging::error(!stream.fail(),
                "IMPERFECTION: malformed RES loadcase in ", filename);

            frame_count    = 0;
            selected_frame = 0;
            read_field     = false;
            numbered_frames.clear();
            processed_frames.clear();
            continue;
        }

        // Resolve frame identity from a NODE displacement field header
        if (line.rfind("FIELD,", 0) == 0) {
            read_field     = false;
            selected_frame = 0;

            if (current_step != step) continue;

            const std::size_t name_begin = line.find("NAME=");
            const std::size_t name_end   = line.find(',', name_begin);
            const bool        nodal      = line.find("TYPE=NODE,") != std::string::npos;

            if (!nodal || name_begin == std::string::npos || name_end == std::string::npos) {
                continue;
            }

            const std::string name = line.substr(name_begin + 5, name_end - name_begin - 5);
            std::string       base = name;
            std::string       suffix;

            // Separate an optional numeric frame suffix from the semantic field name
            const std::size_t separator = name.find_last_of('_');
            if (separator != std::string::npos && separator + 1 < name.size()
                && std::all_of(name.begin() + separator + 1, name.end(),
                               [](unsigned char ch) { return std::isdigit(ch) != 0; })) {
                base   = name.substr(0, separator);
                suffix = name.substr(separator + 1);
            }

            const bool displacement = (base == "DISPLACEMENT")
                                   || (base == "BUCKLING_MODE")
                                   || (base == "MODE_SHAPE");

            if (!displacement) continue;

            // The suffix groups aliases belonging to the same written frame
            if (suffix.empty()) {
                selected_frame = ++frame_count;
            } else {
                const auto [it, inserted] = numbered_frames.emplace(suffix, frame_count + 1);
                if (inserted) ++frame_count;
                selected_frame = it->second;
            }

            // Read at most one displacement representation of each selected frame
            read_field = scales.find(selected_frame) != scales.end()
                      && processed_frames.insert(selected_frame).second;
            continue;
        }

        // Every indexed RES field ends with an explicit field terminator
        if (line == "END FIELD") {
            read_field = false;
            continue;
        }

        if (!read_field) continue;

        // Parse the first three Cartesian displacement components of a nodal row
        std::istringstream stream(line);
        std::string        node;
        Vec3               displacement = Vec3::Zero();

        stream >> node >> displacement(0) >> displacement(1) >> displacement(2);
        const bool valid_row = !stream.fail() && displacement.allFinite();
        const ID   id        = valid_row ? model.compiled_node_id(node) : -1;

        // Validate the entire row before accumulating the geometric contribution
        logging::error(valid_row,
            "IMPERFECTION: malformed RES displacement in ", filename);
        logging::error(id >= 0 && static_cast<std::size_t>(id) < offsets.size(),
            "IMPERFECTION: unknown RES node ", node);

        offsets[static_cast<std::size_t>(id)] += scales.at(selected_frame) * displacement;
        found_frames.insert(selected_frame);
    }

    // Ensure every requested frame has supplied at least one valid nodal value
    for (const auto& [frame, scale] : scales) {
        (void) scale;
        logging::error(found_frames.find(frame) != found_frames.end(),
            "IMPERFECTION: RES displacement frame ", frame, " not found in STEP=", step);
    }
}

/**
 * Reads selected nodal displacement frames from an ASCII FEMaster FRD file.
 *
 * FRD 1PSTEP metadata contains the explicit analysis-step and frame numbers,
 * so unlike RES no field-name-based reconstruction is required. Only the first
 * nodal DISP record for each selected step/frame pair is used. The FRD node
 * numbering scheme encodes its instance and part-local node identifiers as
 * 100000000 * instance_id + local_id.
 *
 * @param filename Source FRD path.
 * @param step Analysis step to read.
 * @param scales Requested one-based frame numbers and scaling factors.
 * @param model Current compiled model providing the FRD node numbering.
 * @param offsets Accumulated global Cartesian nodal imperfections.
 */
void read_frd(const std::string&           filename,
              int                          step,
              const std::map<int, Precision>& scales,
              model::Model&                model,
              std::vector<Vec3>&           offsets) {
    std::ifstream file(filename);

    // Validate the requested result source
    logging::error(file.is_open(),
        "IMPERFECTION: cannot open FRD file ", filename);

    // Reconstruct the same external node labels used by the FRD writer
    std::unordered_map<ID, ID> nodes;
    for (std::size_t row = 0; row < model._data->node_mapping.size(); ++row) {
        const auto& [instance, local] = model._data->node_mapping[row];

        logging::error(instance != nullptr,
            "IMPERFECTION: invalid compiled node mapping");

        const ID external = static_cast<ID>(100000000) * instance->instance_id + local;
        nodes.emplace(external, static_cast<ID>(row));
    }

    // Track active FRD block metadata and avoid importing duplicate DISP fields
    std::string   line;
    std::set<int> processed_frames;
    std::set<int> found_frames;

    int  current_step  = -1;
    int  current_frame = -1;
    bool read_field    = false;

    while (std::getline(file, line)) {
        // A 1PSTEP record identifies the frame of the following result block
        const std::size_t step_marker = line.find("1PSTEP");
        if (step_marker != std::string::npos) {
            std::istringstream stream(line.substr(step_marker + 6));
            int block = 0;
            stream >> block >> current_frame >> current_step;

            logging::error(!stream.fail() && current_step > 0 && current_frame > 0,
                "IMPERFECTION: malformed FRD 1PSTEP record in ", filename);

            read_field = false;
            continue;
        }

        // -4 identifies the physical field and the FRD node/element data type
        if (line.rfind(" -4", 0) == 0) {
            std::istringstream stream(line.substr(3));
            std::string        name;
            int                components = 0;
            int                data_type  = 0;

            stream >> name >> components >> data_type;
            logging::error(!stream.fail(),
                "IMPERFECTION: malformed FRD field header in ", filename);

            const bool displacement = name == "DISP" && components >= 3 && data_type == 1;
            const bool requested    = current_step == step
                                   && scales.find(current_frame) != scales.end();

            read_field = displacement && requested
                      && processed_frames.insert(current_frame).second;
            continue;
        }

        // -3 terminates the current result field
        if (line.rfind(" -3", 0) == 0) {
            read_field = false;
            continue;
        }

        if (!read_field || line.rfind(" -1", 0) != 0) continue;

        // Parse the first three Cartesian displacement components of a nodal row
        const bool complete = line.size() >= 49;
        const ID   external = complete ? static_cast<ID>(std::stoll(line.substr(3, 10))) : -1;
        const auto node     = nodes.find(external);

        Vec3 displacement = Vec3::Zero();
        if (complete) {
            for (int d = 0; d < 3; ++d) {
                displacement(d) = static_cast<Precision>(std::stod(line.substr(13 + 12 * d, 12)));
            }
        }

        // Reject incomplete or incompatible nodal results before accumulation
        logging::error(complete,
            "IMPERFECTION: truncated FRD displacement record in ", filename);
        logging::error(node != nodes.end(),
            "IMPERFECTION: unknown FRD node ", external);
        logging::error(displacement.allFinite(),
            "IMPERFECTION: non-finite FRD displacement at node ", external);

        offsets[static_cast<std::size_t>(node->second)] += scales.at(current_frame) * displacement;
        found_frames.insert(current_frame);
    }

    // Ensure each requested frame has provided displacement rows
    for (const auto& [frame, scale] : scales) {
        (void) scale;
        logging::error(found_frames.find(frame) != found_frames.end(),
            "IMPERFECTION: FRD displacement frame ", frame, " not found in STEP=", step);
    }
}

} // namespace

/**
 * Registers model-level geometric imperfections as changes to reference positions.
 *
 * The parser executes the keyword after assembly compilation but before any
 * geometry-dependent reference data are constructed. Direct rows contain global
 * Cartesian nodal offsets. FILE rows select one-based written displacement
 * frames and their scaling factors in an existing FEMaster RES or FRD result.
 *
 * Every contribution changes both POSITION_REFERENCE and the initial POSITION,
 * leaving the structural displacement unknown equal to zero. Multiple keyword
 * occurrences accumulate their geometric offsets independently.
 *
 * @param registry Input grammar being registered.
 * @param model Compiled model whose reference coordinates are modified.
 * @param parser Owner of the input filename for relative FILE resolution.
 */
void register_imperfection(fem::io::dsl::Registry& registry, model::Model& model, Parser& parser) {
    registry.command("IMPERFECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Modify reference coordinates by direct offsets or scaled RES/FRD displacement frames.");
        command.keyword(fem::io::dsl::KeywordSpec::make()
            .key("FILE").optional()
            .key("STEP").optional());

        // Keep the requested frame coefficients for the currently executed keyword
        auto scales = std::make_shared<std::map<int, Precision>>();

        command.on_enter([scales](const fem::io::dsl::Keys& keys) {
            scales->clear();

            const int  step     = (keys.has("STEP") ? keys.get<int>("STEP") : 1);
            const bool has_file = keys.has("FILE");

            // Validate the keyword before processing any nodal or frame rows
            logging::error(step > 0,
                "IMPERFECTION: STEP must be positive");
            logging::error(has_file || !keys.has("STEP"),
                "IMPERFECTION: STEP requires FILE");
            logging::error(!has_file || !keys.raw("FILE").empty(),
                "IMPERFECTION: FILE must not be empty");
        });

        // Direct global Cartesian node offsets
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::negate(fem::io::dsl::Condition::key_present("FILE")))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("NODE")
                    .one<Precision>().name("DX")
                    .one<Precision>().name("DY")
                    .one<Precision>().name("DZ"))
                .bind([&model](const std::string& node, Precision dx, Precision dy, Precision dz) {
                    const ID   id    = model.compiled_node_id(node);
                    const Vec3 delta = Vec3(dx, dy, dz);

                    // Validate the complete nodal contribution before changing geometry
                    logging::error(id >= 0 && static_cast<Index>(id) < model._data->positions_reference->rows,
                        "IMPERFECTION: unknown node ", node);
                    logging::error(delta.allFinite(),
                        "IMPERFECTION: non-finite nodal offset");

                    // Apply the same stress-free coordinate change to both position fields
                    for (Dim d = 0; d < 3; ++d) {
                        (*model._data->positions_reference)(id, d) += delta(d);
                        (*model._data->positions)(id, d) += delta(d);
                    }
                })));

        // Frame numbers are one-based indices in the written result step
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_present("FILE"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<int>().name("FRAME")
                    .one<Precision>().name("SCALE"))
                .bind([scales](int frame, Precision scale) {
                    logging::error(frame > 0 && std::isfinite(scale),
                        "IMPERFECTION: invalid frame or scale");

                    (*scales)[frame] += scale;
                })));

        command.on_exit([&model, &parser, scales](const fem::io::dsl::Keys& keys) {
            if (!keys.has("FILE")) return;

            // Resolve the source filename relative to the original input deck
            std::filesystem::path filename(keys.raw("FILE"));
            if (filename.is_relative()) {
                filename = std::filesystem::path(parser.input_path()).parent_path() / filename;
            }

            // The DSL normalizes keyword values to uppercase. Recover the original
            // directory-entry spelling when the filesystem is case-sensitive.
            const auto uppercase = [](std::string value) {
                std::transform(value.begin(), value.end(), value.begin(),
                    [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
                return value;
            };

            if (!std::filesystem::exists(filename)) {
                std::filesystem::path resolved = filename.has_root_path()
                    ? filename.root_path() : std::filesystem::path(".");

                for (const auto& component : filename.relative_path()) {
                    const std::filesystem::path candidate = resolved / component;
                    if (std::filesystem::exists(candidate)) {
                        resolved = candidate;
                        continue;
                    }

                    // Search only the current path component for a unique case-insensitive match
                    logging::error(std::filesystem::is_directory(resolved),
                        "IMPERFECTION: result path not found: ", filename.string());

                    std::filesystem::path match;
                    const std::string target = uppercase(component.string());

                    for (const auto& entry : std::filesystem::directory_iterator(resolved)) {
                        if (uppercase(entry.path().filename().string()) == target) {
                            logging::error(match.empty(),
                                "IMPERFECTION: ambiguous filename ", component.string());
                            match = entry.path();
                        }
                    }

                    logging::error(!match.empty(),
                        "IMPERFECTION: result file not found: ", filename.string());
                    resolved = match;
                }
                filename = resolved;
            }

            // Select the matching native FEMaster result reader
            const std::string extension = uppercase(filename.extension().string());
            const int         step      = (keys.has("STEP") ? keys.get<int>("STEP") : 1);

            logging::error(extension == ".RES" || extension == ".FRD",
                "IMPERFECTION: FILE must be a RES or FRD result");

            std::vector<Vec3> offsets(
                static_cast<std::size_t>(model._data->positions_reference->rows), Vec3::Zero());

            if (extension == ".RES") {
                read_res(filename.string(), step, *scales, model, offsets);
            } else {
                read_frd(filename.string(), step, *scales, model, offsets);
            }

            // Commit the accumulated imperfection as a stress-free reference geometry
            auto& reference = *model._data->positions_reference;
            auto& current   = *model._data->positions;

            for (Index node = 0; node < reference.rows; ++node) {
                for (Dim d = 0; d < 3; ++d) {
                    const Precision delta = offsets[static_cast<std::size_t>(node)](d);
                    reference(node, d) += delta;
                    current(node, d)   += delta;
                }
            }
        });
    });
}

} // namespace fem::io::reader::commands
