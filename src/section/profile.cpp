//
// Created by f_eggers on 22.11.2024.
//

#include "profile.h"
#include "../core/logging.h"

fem::Profile::Profile(const std::string& name, Precision A, Precision I_y, Precision I_z, Precision I_t)
    : Namable(name)
    , A(A)
    , I_y(I_y)
    , I_z(I_z)
    , I_t(I_t) {}

void fem::Profile::info() {
    logging::info(true, "Profile: ", name);
    logging::info(true, "   A  : ", A);
    logging::info(true, "   I_y: ", I_y);
    logging::info(true, "   I_z: ", I_z);
    logging::info(true, "   I_t: ", I_t);
}