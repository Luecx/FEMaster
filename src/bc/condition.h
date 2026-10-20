/**
 * @file condition.h
 * @brief Defines the common interface for FEMaster model conditions.
 *
 * A condition describes externally prescribed behavior in the currently
 * active model state. Structural supports, structural loads and thermal boundary
 * conditions all share this semantic layer even though they contribute to
 * different parts of the discrete finite-element system.
 *
 * Every condition receives the complete set of solver-facing assembly outputs:
 *
 *     f   : nodal right-hand-side field,
 *     C,d : algebraic constraint equations,
 *     K_c : condition-dependent matrix contribution.
 *
 * A concrete condition modifies only the outputs required by its physical
 * definition and leaves the remaining outputs unchanged. This keeps the common
 * interface general without introducing separate Dirichlet, Neumann, Mixed or
 * load-capability hierarchies.
 *
 * Condition history remains independent of assembly. The load case owns independent initial and target condition states, while
 * the reader resolves original identifiers and semantic replacement rules.
 * Named collectors retain reusable definitions independently of the active set.
 *
 * @see ConditionState
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#pragma once

#include "amplitude.h"

#include "../constraints/types/equation.h"
#include "../core/printable.h"
#include "../core/types_eig.h"
#include "../core/types_num.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::model {
struct Field;
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common polymorphic definition of a FEMaster condition.
 *
 * A Condition represents one physical prescription such as a structural
 * support, concentrated load or thermal convection boundary. LoadCase owns two ConditionState collections for the initial and target
 * physical prescriptions. Named
 * collectors retain reusable definitions, while the reader decides which
 * definitions participate in an analysis. The solver does not own additional
 * condition history or look up original input identifiers.
 *
 * The reader replaces active definitions based on their original input targets,
 * prescribed component masks and resolved regions. These keyword-specific
 * comparisons and all original input identifiers belong entirely to the
 * reader; the physical Condition remains independent of input syntax.
 *
 * ConditionState membership is based solely on shared pointer identity;
 * distinct condition objects may intentionally address overlapping regions.
 *
 * `apply()` is the single solver-facing assembly interface. Every call receives
 * all available outputs:
 *
 *     rhs             external structural or thermal source field,
 *     equations       algebraic rows C x = d,
 *     system_dof_ids  active equation numbering used by matrix terms,
 *     matrix          condition-dependent sparse operator or tangent entries.
 *
 * Concrete conditions write only the contributions they own. A structural load
 * modifies `rhs`, a support modifies `equations`, and convection modifies both
 * the thermal `rhs` and `matrix`. Passing all outputs on every call avoids
 * capability-specific base classes and lets future conditions contribute to any
 * combination without changing the common interface.
 *
 * amplitude_ optionally references a reusable scalar time history. Concrete
 * conditions decide which of their physical parameters it scales and how the
 * supplied evaluation time is interpreted. The flag ignore_amplitude allows
 * assembly paths such as spatial load-basis construction to bypass this scaling.
 * Numerical scaling and amplitude evaluation belong to the assembly caller, not
 * to the physical Condition definition.
 *
 * No active flag, parser-owned lookup table, step-transition state or numerical
 * commit/rollback state is stored in this base class. Such concerns belong to
 * the input reader and the appropriate numerical analysis. The concrete
 * condition classes store their own physical magnitudes and geometry.
 */
struct Condition : Printable {
    // Shared ownership retains physical definitions in model history or reusable
    // named groups without coupling their lifetime to an analysis object.
    using Ptr = std::shared_ptr<Condition>;

    // Shared time history evaluated by concrete conditions that support
    // amplitude scaling. A null pointer denotes no explicit amplitude.
    Amplitude::Ptr amplitude_ = nullptr;

    // Polymorphic lifetime shared between active storage and named collectors
    ~Condition() override = default;

    // Assemble the contribution of this condition to the global finite
    // element system at the specified time and step progress.
    //
    // All structural and thermal conditions use this common interface.
    // Depending on the concrete condition type, the implementation may
    // contribute to one or more of the following outputs:
    //
    // - rhs:       Accumulates external nodal forces, moments or thermal
    //              loads in the corresponding global field.
    //
    // - equations: Appends prescribed displacement or temperature equations
    //              of the form C*u = d. These are subsequently handled by
    //              the constraint transformation.
    //
    // - matrix:    Appends triplets representing contributions to the global
    //              system matrix, such as thermal convection terms.
    //
    // The outputs are modified in place. Contributions are accumulated rather
    // than assigned, allowing multiple conditions to act on the same nodes
    // or degrees of freedom. Outputs not relevant to the concrete condition
    // remain unchanged.
    //
    // model_data provides the compiled model geometry, regions and other
    // information required to evaluate the physical condition.
    // system_dof_ids maps nodal degrees of freedom to global system indices
    // when matrix contributions are assembled.
    //
    // Every derived condition is responsible for evaluating its own physical
    // contribution, including coordinate transformations and integration
    // over the relevant nodes, surfaces or elements.
    // Assemble this immutable physical prescription using a scalar contribution
    // factor chosen by the analysis. Supports and prescribed temperatures do not
    // scale their Dirichlet values; load conditions scale the complete resultant.
    // Amplitude sampling and step interpolation belong to the caller.
    virtual void apply(
        model::ModelData&      model_data,
        model::Field&          rhs,
        constraint::Equations& equations,
        const SystemDofIds&    system_dof_ids,
        TripletList&           matrix,
        Precision              scale = Precision(1)
    ) = 0;

    // Diagnostics
    std::string str() const override = 0;
};

} // namespace fem::bc
