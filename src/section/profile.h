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
    // Optional product of inertia around local y-z (about the element axis x).
    // Convention: I_yz = integral_A(y*z*dA), i.e. no leading minus sign.
    Precision I_yz;
    // Offsets relative to SMP (shear/mass reference axis).
    // e = SP - SMP, ref = REF - SMP
    Precision e_y;
    Precision e_z;
    Precision ref_y;
    Precision ref_z;

    // constructor
    Profile(const std::string& name,
            Precision A,
            Precision I_y,
            Precision I_z,
            Precision I_t,
            Precision I_yz = 0,
            Precision e_y = 0,
            Precision e_z = 0,
            Precision ref_y = 0,
            Precision ref_z = 0);

    using Ptr = std::shared_ptr<Profile>;

    void info();
};

}



#endif //PROFILE_H
