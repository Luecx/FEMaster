#pragma once

#include <string>
#include <vector>

#include "equation.h"

namespace fem::constraint {

struct EquationFormatOptions {
    int precision = 3;
    bool show_sign = true;
};

std::string format_equation(const Equation& equation,
                            const EquationFormatOptions& opt = {});
std::vector<std::string> format_equations(const Equations& equations,
                                          const EquationFormatOptions& opt = {});

} // namespace fem::constraint

