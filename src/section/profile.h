//
// Created by f_eggers on 22.11.2024.
//

#ifndef PROFILE_H
#define PROFILE_H

#include "../core/types_eig.h"
#include "../data/namable.h"

namespace fem {

struct Profile : public Namable {

    Precision A;

    Precision I_y;
    Precision I_z;
    Precision I_t;

    // constructor
    Profile(const std::string& name, Precision A, Precision I_y, Precision I_z, Precision I_t);

    using Ptr = std::shared_ptr<Profile>;

    void info();
};

}



#endif //PROFILE_H
