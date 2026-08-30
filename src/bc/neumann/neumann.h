/**
 * @file neumann.h
 * @brief Declares the common interface for natural boundary conditions.
 *
 * Neumann conditions contribute prescribed generalized fluxes to a model field.
 * Structural forces, tractions and inertia loads use the same contract as
 * thermal surface heat flux and the source part of convection. Optional
 * coordinate systems and amplitudes are shared by all concrete conditions.
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
 * @brief Common abstraction for prescribed natural boundary conditions.
 *
 * A Neumann condition converts its physical definition into contributions to a
 * nodal model field. Structural conditions assemble generalized forces into the
 * usual mechanical components. Thermal conditions reuse the existing
 * three-component surface integration and store scalar heat-flow contributions
 * exclusively in component zero.
 *
 * The optional coordinate system defines vector-valued conditions in a local
 * basis and the optional amplitude scales the nominal condition at analysis
 * time. Concrete conditions own their target regions and physical parameters.
 */
struct Neumann : BoundaryCondition, Printable {
    using Ptr = std::shared_ptr<Neumann>;

    cos::CoordinateSystem::Ptr orientation_ = nullptr;
    Amplitude::Ptr amplitude_ = nullptr;

    virtual ~Neumann() = default;

    virtual void apply(model::ModelData& model_data,
                       model::Field&     bc,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;

    std::string str() const override = 0;
};

} // namespace fem::bc
