/**
 * @file cload_common.h
 * @brief Shared concentrated-load data syntax and model construction for both readers.
 */
#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "../../../bc/neumann/load_c.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands::cload_common {

inline constexpr const char* missing_token = "__CLOAD_MISSING__";
inline constexpr const char* empty_token   = "__CLOAD_EMPTY__";

inline Precision magnitude(const std::string& token) {
    if (token == empty_token || token == missing_token) {
        throw std::runtime_error("CLOAD: missing numerical component");
    }
    std::size_t consumed = 0;
    const double value = std::stod(token, &consumed);
    if (consumed != token.size()) {
        throw std::runtime_error("CLOAD: invalid numerical component '" + token + "'");
    }
    return static_cast<Precision>(value);
}

/**
 * A data row has either TARGET,DOF,MAGNITUDE or TARGET,Fx,Fy,Fz[,Mx,My,Mz].
 * The DSL supplies separate sentinels for omitted tail entries and explicit
 * empty fields. This preserves the original data-row length after normalization,
 * allowing the two forms to be mixed within one CLOAD keyword block.
 */
inline Vec6 parse(const std::array<std::string, 6>& tokens) {
    std::size_t count = 0;
    while (count < tokens.size() && tokens[count] != missing_token) ++count;
    for (std::size_t i = count; i < tokens.size(); ++i) {
        if (tokens[i] != missing_token) {
            throw std::runtime_error("CLOAD: values cannot follow omitted tail components");
        }
    }

    Vec6 values = Vec6::Zero();
    if (count == 2) {
        if (tokens[0] == empty_token) {
            throw std::runtime_error("CLOAD: missing DOF");
        }
        std::size_t consumed = 0;
        const int dof = std::stoi(tokens[0], &consumed);
        if (consumed != tokens[0].size() || dof < 1 || dof > 6) {
            throw std::runtime_error("CLOAD: DOF must be an integer in [1,6]");
        }
        values[dof - 1] = magnitude(tokens[1]);
        return values;
    }

    if (count < 3 || count > 6) {
        throw std::runtime_error("CLOAD: expected DOF,magnitude or at least Fx,Fy,Fz");
    }

    for (std::size_t i = 0; i < count; ++i) {
        values[static_cast<Index>(i)] =
            tokens[i] == empty_token ? Precision(0) : magnitude(tokens[i]);
    }
    return values;
}

inline void add(model::Model& model,
                model::NodeRegion::Ptr region,
                const Vec6& values,
                cos::CoordinateSystem::Ptr orientation,
                bc::Amplitude::Ptr amplitude) {
    auto load = std::make_shared<bc::CLoad>();
    load->region_      = std::move(region);
    load->values_      = values;
    load->orientation_ = std::move(orientation);
    load->amplitude_   = std::move(amplitude);
    model.add_load(std::move(load));
}

} // namespace fem::io::reader::commands::cload_common
