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

    std::ostream& operator<<(std::ostream& os) const {
        os << "Profile: ";
        os << "   A  : " << A << " ";
        os << "   I_y: " << I_y << " ";
        os << "   I_z: " << I_z << " ";
        os << "   J_t: " << J_t << " ";
        return os;
    }
};

}



#endif //PROFILE_H
