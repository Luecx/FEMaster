/**
 * @file condition_manager.h
 * @brief Declares the model-owned collection of currently active conditions.
 *
 * Conditions describe physical contributions to the finite-element equations:
 * structural loads contribute to the external force vector, prescribed
 * displacements contribute constraint equations, and thermal conditions may
 * contribute source terms, boundary operators or temperature constraints.
 * Their numerical behavior is implemented by the individual Condition classes,
 * not by ConditionManager.
 *
 * ConditionManager only records which physical definitions are currently
 * available for assembly. It stores shared condition pointers in independent
 * input-history families. A family represents the scope of an input operation
 * such as CLOAD, DLOAD, BOUNDARY or CONVECTION; it is not necessarily a unique
 * concrete C++ condition type. Keeping these families separate permits the
 * reader to reset one input-history scope without affecting unrelated loads
 * or boundary conditions.
 *
 * Identification and replacement belong exclusively to the input reader.
 * The reader resolves original input identifiers, handles OP=MOD and OP=NEW,
 * and compares targets using reader-side region and DOF-mask rules. In particular,
 * one input definition may produce several conditions after nodal TRANSFORM
 * processing. The manager does not store identifiers, resolve matching targets
 * or infer which physical definitions should replace each other.
 *
 * Membership is determined by the identity of the shared pointer, not by its
 * numerical contents or the region it addresses. Consequently two different
 * objects may contribute to the same target, whereas inserting the same pointer
 * twice into one family is idempotent. Pointers remain valid while referenced by
 * the manager, a reusable named collector or the input reader. Removing a
 * pointer from the manager only removes that manager-owned reference.
 *
 * The manager contains neither amplitude interpolation nor step-time,
 * displacement, Newton or material state. Model assembly traverses the active
 * families and invokes Condition::apply() on their stored definitions.
 *
 * @see Condition
 * @see LoadCollector
 * @see SupportCollector
 * @see ThermalCollector
 *
 * @author Finn Eggers
 * @date 08.10.2026
 */

#pragma once

#include "condition.h"

#include <array>
#include <cstddef>
#include <unordered_set>

namespace fem::bc {

/**
 * @enum ConditionFamily
 * @brief Identifies the independent input-history scope of a condition.
 *
 * A family defines which stored conditions are affected by a reader-side
 * replacement or reset operation. It does not prescribe the assembly operator:
 * for example, DLOAD and DSLOAD may each produce different concrete structural
 * load classes, and physically similar loads may belong to different families
 * because their input commands have independent history.
 *
 * Model assembly selects the relevant families explicitly when constructing
 * structural forces, structural supports and thermal contributions. No condition
 * is classified through its dynamic C++ type by ConditionManager.
 *
 * N_CONDITION_FAMILIES is the array extent, not a valid condition family.
 */
enum ConditionFamily : std::size_t {
    // Prescribed structural translations and rotations, including BOUNDARY
    SUPPORT = 0,

    // Concentrated nodal forces and moments
    CLOAD,

    // Element-addressed distributed loads
    DLOAD,

    // Surface-addressed distributed loads
    DSLOAD,

    // FEMaster-native scalar surface pressures
    PLOAD,

    // FEMaster-native distributed volume forces
    VLOAD,

    // Prescribed rigid-body accelerations and inertia loads
    INERTIAL_LOAD,

    // Prescribed thermal degrees of freedom
    TEMPERATURE,

    // Applied surface heat-flow densities
    HEAT_FLUX,

    // Convective thermal boundary conditions
    CONVECTION,

    // Array extent; always keep this final
    N_CONDITION_FAMILIES
};

/**
 * @class ConditionManager
 * @brief Owns references to the active conditions of one finite-element model.
 *
 * Each family is an independent unordered set of shared Condition pointers.
 * The outer array uses ConditionFamily as its index, so removing or clearing
 * entries in one family leaves every other family unchanged.
 *
 * A set is appropriate because assembly requires iteration over all active
 * physical definitions, while the input reader must also check membership and
 * remove a previously resolved pointer efficiently. The set uses the shared
 * pointer's stored address for hashing and equality, rather than evaluating
 * semantic input targets. Distinct definitions on overlapping regions are
 * therefore preserved and their numerical contributions remain additive.
 *
 * The manager does not own a separate identifier index or assign numeric
 * condition IDs. The reader keeps its own mapping from original identifiers
 * to the corresponding shared pointers, including multiple physical fragments
 * created from one input target. The reader must maintain that index when it
 * removes or clears active definitions; this manager intentionally cannot
 * synchronize external lookup tables.
 *
 * Named collectors retain reusable definitions independently of the active
 * state. The reader may insert their pointers for an analysis and remove those
 * newly inserted references afterwards. The Boolean result of add() allows
 * the reader to distinguish a new insertion from an already active pointer.
 *
 * No iteration order is guaranteed by the unordered sets. Assembly results
 * may therefore differ by floating-point rounding between iteration orders,
 * but the mathematical set of active contributions remains the same.
 *
 * No step history, interpolation, semantic replacement or solver state is
 * committed here; those responsibilities remain with the reader, the concrete
 * condition definitions and the numerical analysis.
 */
class ConditionManager {
public:
    // Physical definition pointers, unique by pointer identity within a family
    using Conditions = std::unordered_set<Condition::Ptr>;

private:
    // Independently active definitions for every input-history family
    std::array<Conditions, N_CONDITION_FAMILIES> conditions_;

public:
    // Insert a physical condition reference. Repeated insertion of the same
    // pointer into one family has no effect; other families remain independent.
    bool add(ConditionFamily family, Condition::Ptr condition);

    // Remove or test exactly the specified pointer within one family. Neither
    // operation compares regions, prescribed components or numerical values.
    bool remove(ConditionFamily family, const Condition::Ptr& condition);
    bool contains(ConditionFamily family, const Condition::Ptr& condition) const;

    // Read the live active set without copying it. Readers must not retain
    // iterators across later insertions or removals of this family.
    const Conditions& get(ConditionFamily family) const;
};

} // namespace fem::bc
