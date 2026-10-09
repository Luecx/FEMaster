/**
 * @file c3d5.h
 * @brief Declares the five-node pyramid as a degenerate C3D8 solid.
 *
 * The four pyramid base nodes occupy the lower hexahedral face. All four
 * upper C3D8 corners reference the same fifth (apex) node. The inherited
 * solid formulation therefore assembles the pyramid through eight local
 * connectivity slots without introducing a separate element formulation.
 *
 * @see C3D8
 * @see SolidElement
 */

#pragma once

#include "c3d8.h"

namespace fem::model {

/**
 * @brief Five-node pyramid represented by a collapsed eight-node hexahedron.
 *
 * Constructor connectivity is [base 1, base 2, base 3, base 4, apex].
 * The expanded eight-node connectivity remains owned by C3D8.
 */
struct C3D5 : C3D8 {
    C3D5(ID elem_id, const std::array<ID, 5>& node_ids);

    // Reconstruct from the five independent connectivity entries so copying
    // during Model::compile() preserves the concrete pyramid element type.
    ElementPtr copy() const override;
    std::string type_name() const override { return "C3D5"; }
};

} // namespace fem::model
