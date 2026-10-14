/**
 * @file c3d5.cpp
 * @brief Implements the degenerate-C3D8 five-node pyramid connectivity.
 *
 * The pyramid base maps to the lower C3D8 face and the apex is repeated in
 * every upper-corner slot. Shape functions, integration and mechanics are
 * inherited without modification.
 *
 * @see C3D5
 * @see C3D8
 */

#include "c3d5.h"

#include "../geometry/surface/surface3.h"
#include "../geometry/surface/surface4.h"

namespace fem::model {

/**
 * @brief Expands the five pyramid nodes into collapsed hexahedral connectivity.
 *
 * @param elem_id Element identifier.
 * @param node_ids Four base corners followed by the pyramid apex.
 */
C3D5::C3D5(ID elem_id, const std::array<ID, 5>& node_ids)
    : C3D8(elem_id, {
        node_ids[0], node_ids[1], node_ids[2], node_ids[3],
        node_ids[4], node_ids[4], node_ids[4], node_ids[4]
    }) {}

/**
 * Extracts the five physical pyramid faces using Abaqus S1-S5 numbering.
 *
 * The inherited C3D8 connectivity repeats the apex in four upper slots.
 * Returning its S2 face would create a zero-area quadrilateral, so the
 * pyramid instead exposes its base and four properly oriented triangles.
 */
SurfacePtr C3D5::surface(ID surface_id) {
    switch (surface_id) {
        case 1: return std::make_shared<Surface4>(std::array<ID, 4>{node_ids[0], node_ids[3], node_ids[2], node_ids[1]});
        case 2: return std::make_shared<Surface3>(std::array<ID, 3>{node_ids[0], node_ids[1], node_ids[4]});
        case 3: return std::make_shared<Surface3>(std::array<ID, 3>{node_ids[1], node_ids[2], node_ids[4]});
        case 4: return std::make_shared<Surface3>(std::array<ID, 3>{node_ids[2], node_ids[3], node_ids[4]});
        case 5: return std::make_shared<Surface3>(std::array<ID, 3>{node_ids[3], node_ids[0], node_ids[4]});
        default: return nullptr;
    }
}

/**
 * @brief Copies the persistent five-node topology as another C3D5.
 *
 * The expanded connectivity is reconstructed by the constructor. Compilation
 * subsequently remaps all eight references to their dense assembly node ids.
 */
ElementPtr C3D5::copy() const {
    return std::make_shared<C3D5>(elem_id, std::array<ID, 5> {
        node_ids[0], node_ids[1], node_ids[2], node_ids[3], node_ids[4]
    });
}

} // namespace fem::model
