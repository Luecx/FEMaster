/**
 * @file part.h
 * @brief Defines reusable finite-element part data before model compilation.
 *
 * A `Part` owns local topology, named regions and section assignments that may
 * be instantiated multiple times in one `Model`. Node and element identifiers
 * are local to the part and therefore need not form dense ranges. Surfaces and
 * lines retain the same local node connectivity, while sections reference local
 * element regions.
 *
 * `Instance` adds only placement information and never duplicates the part
 * definition. `Model::compile()` expands every instance into the dense,
 * solver-facing `ModelData` representation used by elements and load cases.
 *
 * @see Instance
 * @see Model
 * @see ModelData
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../data/namable.h"
#include "../data/region.h"
#include "../data/sets.h"
#include "../section/section.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fem::model {

/**
 * @brief Reusable local finite-element definition referenced by model instances.
 *
 * A part stores nodes, elements, geometric boundary entities, named regions and
 * section assignments in its own identifier space. These identifiers remain
 * stable while the model is being constructed and may be sparse or otherwise
 * unrelated to the dense solver indices generated later.
 *
 * Sections stored by the part reference local element regions. During
 * compilation each instance receives independent compiled assignments with the
 * affected element identifiers remapped into assembly space.
 *
 * The class intentionally contains no compilation or solver logic. Expansion,
 * coordinate transformation and dense enumeration are responsibilities of
 * `Model::compile()`.
 */
struct Part : public Namable {
    // Shared ownership type used by the model dictionary and by every instance
    // that embeds this part into the assembled model.
    using Ptr = std::shared_ptr<Part>;

    // Part-local topology. Node coordinates are expressed in the local part
    // coordinate system, and all element, surface and line connectivity refers
    // to identifiers from these local containers. Model compilation copies and
    // rewires the entities into dense assembly storage without modifying this
    // reusable definition.
    std::unordered_map<ID, Vec3>       nodes;
    std::unordered_map<ID, ElementPtr> elements;
    std::unordered_map<ID, SurfacePtr> surfaces;
    std::unordered_map<ID, LinePtr>    lines;

    // Named regions in the same local identifier space as the topology above.
    // Each registry also owns its aggregate *ALL region. Compilation creates an
    // instance-qualified assembly region for every local region and maps all
    // contained identifiers through the corresponding instance.
    Sets<NodeRegion>    node_sets   {SET_NODE_ALL};
    Sets<ElementRegion> elem_sets   {SET_ELEM_ALL};
    Sets<SurfaceRegion> surface_sets{SET_SURF_ALL};
    Sets<LineRegion>    line_sets   {SET_LINE_ALL};

    // Section assignments defined on local element regions. Every instance
    // receives an independent compiled section whose region and any spatial
    // properties are transformed or remapped into assembly space.
    std::vector<Section::Ptr> sections;

    // Construction of an empty named part. Topology, regions and sections are
    // populated by the parser or public model interface before the one-way model
    // compilation pass.
    explicit Part(std::string name)
        : Namable(std::move(name)) {}
};

} // namespace fem::model
