/**
 * @file InputDeckParser.cpp
 * @brief Implementation of InputDeckParser (analyse + process).
 */

#include "InputDeckParser.h"

#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <utility>

// Use the parser library without changing it
using namespace fem::reader2;

namespace fem::reader2 {

// ============================================================================
// ctor / top-level
// ============================================================================
InputDeckParser::InputDeckParser(const std::string& file_path, const std::string& output_path)
    : m_file_path(file_path)
    , m_output_path(output_path)
    , m_writer(output_path) { }

void InputDeckParser::read() {
    // Phase 1: cheap scan
    analyse();

    // Create model from discovered maxima
    m_model = std::make_shared<fem::model::Model>(m_data.highest_node_id     + 1,
                                                  m_data.highest_element_id  + 1,
                                                  m_data.highest_surface_id  + 1);

    // Phase 2: full parse through the new engine
    process();

    // Finalize writer
    m_writer.close();
}

// ============================================================================
// analyse(): simple tolerant scan for maxima (no complex parsing)
// ============================================================================
void InputDeckParser::analyse() {
    File file(m_file_path);
    Line& L = file.next();

    auto to_int = [](const std::string& s)->int {
        try { return std::stoi(s); } catch (...) { return -1; }
    };

    while (L.type() != LineType::EOF_MARK) {
        // advance until keyword or EOF
        while (L.type() != LineType::KEYWORD && L.type() != LineType::EOF_MARK) {
            L = file.next();
        }
        if (L.type() == LineType::EOF_MARK) break;

        const std::string cmd = L.command(); // uppercased by Line

        if (cmd == "NODE") {
            // consume node data lines; first token is node id
            L = file.next();
            while (L.type() == LineType::DATA) {
                const auto& vals = L.values();
                if (!vals.empty()) {
                    int nid = to_int(vals[0]);
                    if (nid >= 0) m_data.highest_node_id = std::max(m_data.highest_node_id, nid);
                }
                L = file.next();
            }
            continue; // L is now non-DATA; loop continues
        }
        if (cmd == "ELEMENT") {
            // consume element data lines; first token is element id
            L = file.next();
            while (L.type() == LineType::DATA) {
                const auto& vals = L.values();
                if (!vals.empty()) {
                    int eid = to_int(vals[0]);
                    if (eid >= 0) m_data.highest_element_id = std::max(m_data.highest_element_id, eid);
                }
                L = file.next();
            }
            continue;
        }
        if (cmd == "SURFACE") {
            // consume surface data lines; first token (if numeric) is surface id
            L = file.next();
            while (L.type() == LineType::DATA) {
                const auto& vals = L.values();
                if (!vals.empty()) {
                    int sid = to_int(vals[0]);
                    if (sid >= 0) m_data.highest_surface_id = std::max(m_data.highest_surface_id, sid);
                }
                L = file.next();
            }
            continue;
        }

        // unknown keyword: step forward
        L = file.next();
    }
}