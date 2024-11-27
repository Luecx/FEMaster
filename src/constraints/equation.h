//
// Created by f_eggers on 25.11.2024.
//

#ifndef EQUATION_H
#define EQUATION_H

#include "../core/core.h"

namespace fem::constraint {

struct EquationEntry {
    ID node_id;
    Dim dof;
    Precision value;
};

struct Equation {
    std::vector<EquationEntry> entries;

    Equation(std::vector<EquationEntry> entries) : entries(std::move(entries)) {};
    Equation(std::initializer_list<EquationEntry> entries) : entries(entries) {};

    using Equations = std::vector<Equation>;

    static TripletList get_triplets(Equations& equations, SystemDofIds& dofs, Index start_row) {
        TripletList triplets;
        Index row = start_row;
        for (auto& eq : equations) {

            // ignore empty equations
            if (eq.entries.empty()) {
                continue;
            }

            // check if any entry is found, if not, no need to increment the row
            bool has_entry = false;
            for (auto& entry : eq.entries) {
                int dof_id = dofs(entry.node_id, entry.dof);
                // check if the DOF is active
                if (dof_id < 0)
                    continue;

                has_entry = true;
                triplets.emplace_back(row, dof_id, entry.value);
            }
            if (has_entry)
                row++;
        }
        return triplets;
    }
};

using Equations = std::vector<Equation>;
using EquationEntries = std::vector<EquationEntry>;

}





#endif //EQUATION_H
