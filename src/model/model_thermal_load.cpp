/**
 * @file model_thermal_load.cpp
 * @brief Implements model-level assembly of thermal boundary loads.
 *
 * Thermal natural boundary conditions contribute to a scalar nodal right-hand
 * side and may additionally append sparse operator terms. Prescribed heat flux
 * contributes only to the right-hand side, while convection contributes both
 * the ambient source and the boundary conductivity matrix.
 *
 * @see Model::build_thermal_load_matrix
 * @see bc::HeatFlux
 * @see bc::Convection
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "model.h"

#include "../bc/neumann/convection.h"
#include "../bc/neumann/heat_flux.h"
#include "../core/logging.h"

#include <memory>
#include <string>

namespace fem::model {

/**
 * Assembles selected thermal load collectors into nodal RHS storage and sparse
 * boundary-operator triplets.
 *
 * The supplied system mapping must contain exactly one scalar temperature
 * component per node. Every selected collector is restricted to thermal natural
 * conditions so mechanical loads cannot be interpreted as scalar heat input.
 * Heat-flux conditions write only to the returned nodal field. Convection also
 * appends its temperature-dependent boundary matrix to `lhs` using the active
 * thermal system identifiers.
 *
 * @param load_sets Names of load collectors included in the thermal system.
 * @param system_dof_ids Node-by-one thermal system mapping.
 * @param lhs Sparse triplet list receiving Robin boundary-operator terms.
 * @param time Evaluation time supplied to load amplitudes.
 * @return One-component nodal thermal right-hand-side field.
 */
Field Model::build_thermal_load_matrix(
    const std::vector<std::string>& load_sets,
    const SystemDofIds&             system_dof_ids,
    TripletList&                    lhs,
    Precision                       time
) {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");
    logging::error(system_dof_ids.rows() == static_cast<Eigen::Index>(_data->positions->rows),
        "Model: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "Model: thermal DOF map must contain exactly one component per node");

    Field rhs{
        "THERMAL_LOAD",
        FieldDomain::NODE,
        _data->field_rows(FieldDomain::NODE),
        1
    };
    rhs.set_zero();

    for (const std::string& name : load_sets) {
        logging::error(_data->load_cols.has(name),
            "Model: thermal load collector ", name, " does not exist");

        const auto collector = _data->load_cols.get(name);
        logging::error(collector != nullptr,
            "Model: thermal load collector ", name, " is not initialized");

        for (const auto& load : collector->entries()) {
            if (!load) continue;

            const bool thermal = std::dynamic_pointer_cast<bc::HeatFlux>(load) != nullptr
                              || std::dynamic_pointer_cast<bc::Convection>(load) != nullptr;
            logging::error(thermal,
                "Model: thermal load collector ", name,
                " contains non-thermal condition ", load->str());

            load->apply(
                *_data,
                rhs,
                time,
                false,
                &system_dof_ids,
                &lhs
            );
        }
    }

    return rhs;
}

} // namespace fem::model
