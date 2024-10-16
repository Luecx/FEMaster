//
// Created by Luecx on 04.09.2023.
//
#include "line.h"

#include "../core/logging.h"

namespace fem::reader {

Line& fem::reader::Line::operator=(const std::string& line) {
    // reset fields
    m_line = line;
    m_values.clear();
    m_keys.clear();
    m_command.clear();
    m_type = END_OF_FILE;

    // Remove leading whitespaces and special characters
    m_line.erase(m_line.begin(), std::find_if(m_line.begin(), m_line.end(), [](int ch) { return !std::isspace(ch); }));

    // Remove trailing whitespaces and special characters
    m_line.erase(std::find_if(m_line.rbegin(), m_line.rend(), [](int ch) { return !std::isspace(ch); }).base(),
                 m_line.end());

    // remove every whitespace
    m_line.erase(std::remove(m_line.begin(), m_line.end(), ' '), m_line.end());

    // check if the m_line is a comment --> starts with //, #, **, !
    if (m_line.substr(0, 2) == "//" || m_line[0] == '#' || m_line.substr(0, 2) == "**" || m_line[0] == '!') {
        m_type = COMMENT;
    } else if (m_line.empty()) {
        m_type = EMPTY_LINE;
    } else {
        // check if it is a keyword line
        if (m_line[0] == '*') {
            m_type = KEYWORD_LINE;
            m_line.erase(std::remove(m_line.begin(), m_line.end(), ' '), m_line.end());
            m_line.erase(0, 1);                                                         // remove the leading '*'
            std::transform(m_line.begin(), m_line.end(), m_line.begin(), ::toupper);    // convert to upper case
            std::stringstream ss(m_line);
            std::string       item;
            if (std::getline(ss, m_command, ',')) {    // get the command name
                // Remove trailing special characters from the command
                m_command.erase(
                    std::find_if(m_command.rbegin(), m_command.rend(), [](int ch) { return !std::isspace(ch); }).base(),
                    m_command.end());

                while (std::getline(ss, item, ',')) {
                    size_t pos = item.find("=");
                    if (pos != std::string::npos) {
                        m_keys[item.substr(0, pos)] = item.substr(pos + 1);
                    } else {
                        m_keys[item] = "";
                    }
                }
            }
        } else {
            m_type = DATA_LINE;
            m_line = std::regex_replace(m_line, std::regex(" +"), "");    // replace multiple spaces with a single space
            std::transform(m_line.begin(), m_line.end(), m_line.begin(), ::toupper);    // convert to upper case
            std::stringstream ss(m_line);
            std::string       item;
            while (std::getline(ss, item, ',')) {
                m_values.push_back(item);
            }
        }
    }

    return *this;
}
const std::string& fem::reader::Line::line() const {
    return m_line;
}
const std::vector<std::string>& fem::reader::Line::values() const {
    return m_values;
}
bool fem::reader::Line::has_key(const std::string& key) const {
    return m_keys.find(key) != m_keys.end();
}
const std::string& fem::reader::Line::command() const {
    return m_command;
}
LineType fem::reader::Line::type() const {
    return m_type;
}
size_t fem::reader::Line::count_values() const {
    logging::error(m_type == DATA_LINE, "The 'count_values' function can only be used with DATA_LINE type.");
    return m_values.size();
}
void fem::reader::Line::eof() {
    this->m_type = END_OF_FILE;
}
}    // namespace fem::reader
