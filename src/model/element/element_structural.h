/**
 * @file element_structural.h
 * @brief Declares the abstract base class for structural FEM elements.
 *
 * Structural elements extend `ElementInterface` with routines for stiffness,
 * mass, volume, load application, and stress/strain recovery.
 *
 * @see src/model/beam/beam.h
 */

#pragma once

#include "element.h"

#include "../../material/stress.h"
#include "../../material/strain.h"

namespace fem {
namespace model {

/**
 * @struct StructuralElement
 * @brief Extends `ElementInterface` with structural-specific operations.
 */
struct StructuralElement : ElementInterface {
    explicit StructuralElement(ID elem_id)
        : ElementInterface(elem_id) {}

    ~StructuralElement() override = default;

    virtual Precision volume() = 0;
    virtual MapMatrix stiffness(Precision* buffer) = 0;
    virtual MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) = 0;
    virtual MapMatrix mass(Precision* buffer) = 0;

    virtual Stresses stress(Field& displacement, std::vector<Vec3>& rst) = 0;
    virtual Strains strain(Field& displacement, std::vector<Vec3>& rst) = 0;
    virtual std::vector<Vec6> section_forces(Field& displacement) {return std::vector<Vec6>();};

    virtual void compute_stress_strain(Field& ip_stress,
                                       Field& ip_strain,
                                       Field& displacement,
                                       int ip_offset) = 0;

    virtual void compute_stress_strain_nodal(Field& displacement,
                                             Field& stress,
                                             Field& strain) = 0;

    virtual void apply_vload(Field& node_loads, Vec3 load) = 0;
    virtual void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) = 0;
    virtual void compute_compliance(Field& displacement, Field& result) = 0;
    virtual void compute_compliance_angle_derivative(Field& displacement, Field& result) = 0;
};

} // namespace model
} // namespace fem
