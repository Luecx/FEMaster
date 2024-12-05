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
    Precision I_t;

    // constructor
    Profile(const std::string& name, Precision A, Precision I_y, Precision I_z, Precision I_t) :
        Namable(name), A(A), I_y(I_y), I_z(I_z), I_t(I_t) {}

    using Ptr = std::shared_ptr<Profile>;

    void info() {
        logging::info(true, "Profile: ", name);
        logging::info(true, "   A  : ", A);
        logging::info(true, "   I_y: ", I_y);
        logging::info(true, "   I_z: ", I_z);
        logging::info(true, "   I_t: ", I_t);
    }
};

}



#endif //PROFILE_H
