//
// Created by f_eggers on 22.11.2024.
//

#ifndef PROFILE_H
#define PROFILE_H

#include "../core/core.h"
#include "../data/namable.h"

namespace fem {

struct Profile : public Namable {

    Precision A;

    Precision I_y;
    Precision I_z;
    Precision J_t;

    // constructor
    Profile(const std::string& name, Precision A, Precision I_y, Precision I_z, Precision J_t) :
        Namable(name), A(A), I_y(I_y), I_z(I_z), J_t(J_t) {}

    using Ptr = std::shared_ptr<Profile>;

    void info() {
        logging::info(true, "Profile: ", name);
        logging::info(true, "   A  : ", A);
        logging::info(true, "   I_y: ", I_y);
        logging::info(true, "   I_z: ", I_z);
        logging::info(true, "   J_t: ", J_t);
    }
};

}



#endif //PROFILE_H
