#pragma once

#include "../core/types_eig.h"
#include "../core/types_num.h"

#include <memory>
#include <string>
#include <vector>

namespace fem::pretension {

enum class State {
    Open,
    Force,
    Hold
};

struct InterfacePair {
    ID side_a = -1;
    ID side_b = -1;
};

struct PretensionSection {
    using Ptr = std::shared_ptr<PretensionSection>;

    std::string name;

    std::string interface_surface_set_a;
    std::string interface_surface_set_b;

    Precision total_force = 0;

    Vec3 axis_direction = Vec3::UnitX();
    std::vector<InterfacePair> interface_pairs;
    State state = State::Open;

    Precision prescribed_value = 0;
    Precision locked_gap = 0;
    Precision last_solved_gap = 0;
    bool has_solved_gap = false;

    bool prepared = false;
};

}
