/**
 * @file instance.h
 * @brief Defines placement and compiled identifier mappings of one part instance.
 *
 * An `Instance` references one reusable part definition together with a rigid
 * translation and rotation. It does not duplicate nodes, elements, sets,
 * surfaces or sections. These entities remain owned by the referenced `Part`
 * and are expanded once by `Model::compile()` into solver-facing dense topology.
 *
 * Instances are retained through compilation in the semantic dictionary owned
 * by `ModelData`; the concrete type is only included and operated by `Model`.
 *
 * The transformation maps local part coordinates to assembled global
 * coordinates according to
 *
 * \f[
 *     \boldsymbol x = \boldsymbol R\,\boldsymbol X + \boldsymbol t.
 * \f]
 *
 * @see Part
 * @see Model
 * @see ModelData
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../core/types_eig.h"
#include "../data/namable.h"
#include "part.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace fem::model {

/**
 * @brief Named rigid placement of one reusable finite-element part.
 *
 * Translation and rotation are the only geometric instance-specific semantic
 * state. The four local-to-dense maps are populated by the single compile and
 * support parser/API operations that keep addressing local identifiers after the
 * semantic topology has been frozen.
 */
struct Instance : public Namable {
    using Ptr = std::shared_ptr<Instance>;

    Part::Ptr part = nullptr;

    Vec3 translation = Vec3::Zero();
    Mat3 rotation    = Mat3::Identity();

    std::unordered_map<ID, ID> node_ids;
    std::unordered_map<ID, ID> element_ids;
    std::unordered_map<ID, ID> surface_ids;
    std::unordered_map<ID, ID> line_ids;

    Instance(std::string name,
             Part::Ptr   part,
             Vec3        translation = Vec3::Zero(),
             Mat3        rotation    = Mat3::Identity())
        : Namable(std::move(name)),
          part       (std::move(part)),
          translation(translation),
          rotation   (rotation) {}
};

} // namespace fem::model
