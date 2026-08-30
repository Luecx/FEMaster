/**
 * @file neumann.h
 * @brief Declares the common interface for natural/load-side boundary conditions.
 *
 * `Neumann` is the common load-side boundary-condition abstraction used by
 * structural and thermal analyses. Every concrete condition receives the same
 * assembly interface: a nodal right-hand-side field is always available, while
 * a system DOF map and sparse left-hand-side triplet list may additionally be
 * supplied by analyses that allow a boundary condition to modify the operator.
 *
 * Most concrete loads, such as concentrated forces, pressure, body forces and
 * prescribed heat flux, contribute only to the right-hand side and therefore
 * ignore the optional system objects. Linear convection is different: although
 * it is mathematically a Robin condition, it is kept in the common load-side
 * hierarchy and uses the same single `apply()` call to assemble both its
 * ambient-temperature source term and its temperature-dependent boundary matrix.
 *
 * @see LoadCollector
 * @see HeatFlux
 * @see Convection
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "../amplitude.h"
#include "../bc.h"

#include "../../core/printable.h"
#include "../../core/types_cls.h"
#include "../../core/types_eig.h"
#include "../../cos/coordinate_system.h"
#include "../../data/field.h"
#include "../../data/region.h"

#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common abstraction for conditions assembled on the load side.
 *
 * The interface deliberately describes algebraic contributions rather than the
 * strict mathematical classification of the boundary condition. A conventional
 * Neumann load modifies only the right-hand side. A Robin-type load may use the
 * optional DOF map and triplet list to modify the left-hand-side operator in the
 * same call. This keeps load application uniform for the analysis while allowing
 * simple loads to remain simple.
 *
 * Structural load fields conventionally use six components per node. Thermal
 * loads reuse the existing vector integration infrastructure and store their
 * scalar heat-flow contribution in component zero of a three-component nodal
 * field; components one and two remain zero.
 */
struct Neumann : BoundaryCondition, Printable {
    using Ptr = std::shared_ptr<Neumann>;

    // Optional coordinate system for vector-valued natural conditions.
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Optional scalar time history applied to the nominal condition.
    Amplitude::Ptr amplitude_ = nullptr;

    virtual ~Neumann() = default;

    /**
     * @brief Assembles this condition into the global load-side system.
     *
     * `rhs` is the mandatory nodal right-hand-side field. Pure Neumann loads
     * modify only this field. Conditions that also contribute to the global
     * operator receive `system_dof_ids` and `lhs` from the calling analysis and
     * append sparse matrix entries to `lhs`. Such conditions must validate that
     * both optional pointers are non-null before using them.
     *
     * The optional objects are pointers rather than references so existing
     * structural assembly paths do not need to manufacture an unused matrix or
     * DOF map for loads that cannot modify the left-hand side.
     *
     * @param model_data Compiled model topology, geometry and material data.
     * @param rhs Nodal right-hand-side field receiving load contributions.
     * @param time Analysis time used to evaluate an attached amplitude.
     * @param ignore_amplitude If true, assemble the nominal spatial load without
     *                         multiplying by the attached amplitude.
     * @param system_dof_ids Optional node-by-component global system DOF map.
     *                       Required only by conditions that assemble `lhs`.
     * @param lhs Optional sparse triplet list receiving left-hand-side entries.
     *            Required only by conditions that modify the system operator.
     */
    virtual void apply(model::ModelData&       model_data,
                       model::Field&           rhs,
                       Precision               time,
                       bool                    ignore_amplitude = false,
                       const SystemDofIds*      system_dof_ids = nullptr,
                       TripletList*             lhs = nullptr) = 0;

    // Return the concrete condition in a compact diagnostic representation.
    std::string str() const override = 0;
};

} // namespace fem::bc
