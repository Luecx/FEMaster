//
// Created by f_eggers on 22.11.2024.
//

#include "profile.h"
#include "../core/logging.h"

fem::Profile::Profile(const std::string& name,
                      Precision A,
                      Precision I_y,
                      Precision I_z,
                      Precision I_t,
                      Precision I_yz,
                      Precision e_y,
                      Precision e_z,
                      Precision ref_y,
                      Precision ref_z)
    : Namable(name)
    , A(A)
    , I_y(I_y)
    , I_z(I_z)
    , I_t(I_t)
    , I_yz(I_yz)
    , e_y(e_y)
    , e_z(e_z)
    , ref_y(ref_y)
    , ref_z(ref_z) {}

void fem::Profile::info() {
    logging::info(true, "Profile: ", name);
    logging::info(true, "   A  : ", A);
    logging::info(true, "   I_y: ", I_y);
    logging::info(true, "   I_z: ", I_z);
    logging::info(true, "   I_t: ", I_t);
    logging::info(true, "   I_yz: ", I_yz);
    logging::info(true, "   e_y: ", e_y);
    logging::info(true, "   e_z: ", e_z);
    logging::info(true, "   ref_y: ", ref_y);
    logging::info(true, "   ref_z: ", ref_z);
}
