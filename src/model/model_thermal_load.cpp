/**
 * @file model_thermal_load.cpp
 * @brief Implements thermal boundary-load assembly.
 */

#include "model.h"

#include "../bc/neumann/neumann.h"
#include "../bc/robin/robin.h"
#include "../core/logging.h"

#include <memory>
#include <string>

namespace fem::model {

/**
 * @brief Assembles selected thermal boundary conditions.
 *
 * Thermal Neumann conditions contribute only to the scalar nodal RHS. Robin
 * conditions contribute to the same RHS and append sparse operator triplets
 * directly through the supplied thermal system DOF mapping.
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

            if (auto neumann = std::dynamic_pointer_cast<bc::ThermalNeumann>(load)) {
                neumann->apply(*_data, rhs, time, false);
                continue;
            }

            if (auto robin = std::dynamic_pointer_cast<bc::Robin>(load)) {
                robin->apply(*_data, rhs, system_dof_ids, lhs, time, false);
                continue;
            }

            logging::error(false,
                "Model: thermal load collector ", name,
                " contains non-thermal condition ", load->str());
        }
    }

    return rhs;
}

} // namespace fem::model
