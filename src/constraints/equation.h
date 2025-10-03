#pragma once
//
// Created by f_eggers on 25.11.2024 (revised: Constraint Equations with RHS)
//

#include <cstdint>
#include <vector>

#include "../core/types_eig.h"
#include "../core/types_cls.h"
#include "../core/types_num.h"

namespace fem::constraint {

    /// One term c * u(node_id, dof)
    struct EquationEntry {
        ID        node_id;   ///< global node id
        Dim       dof;       ///< component index (0..ndim-1)
        Precision coeff;     ///< coefficient c for this dof
    };

    enum class EquationSourceKind : uint8_t {
        Unknown = 0,
        Support,
        Connector,
        Coupling,
        Tie,
        Manual
    };

    /// Linear equation: sum_i coeff_i * u(node_i,dof_i) = rhs
    struct Equation {
        std::vector<EquationEntry> entries;
        Precision rhs = Precision(0);
        EquationSourceKind source = EquationSourceKind::Unknown;
        Index source_index = 0;

        Equation() = default;
        Equation(std::vector<EquationEntry> e, Precision r = 0)
            : entries(std::move(e)), rhs(r) {}
        Equation(std::initializer_list<EquationEntry> e, Precision r = 0)
            : entries(e), rhs(r) {}
    };

    using Equations = std::vector<Equation>;

} // namespace fem::constraint
