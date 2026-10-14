#pragma once

#include "../core/types_eig.h"
#include "../core/types_num.h"

#include <memory>
#include <string>
#include <vector>

namespace fem::pretension {

enum class Control {
    Force,
    Displacement
};

enum class State {
    Open,
    Loading,
    Locked
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

    bool nodal_control = false; // Abaqus conditions own load/boundary history.
    Vec3 normal_override = Vec3::Zero();
    ID control_node = -1; // Dense node ID; only component 0 is active.
    Precision total_force = 0;

    Vec3 axis_direction = Vec3::UnitX();
    std::vector<InterfacePair> interface_pairs;
    Control control = Control::Displacement;
    State state = State::Open;

    Precision prescribed_value = 0;
    Precision locked_gap = 0;
    Precision last_solved_gap = 0;
    bool has_solved_gap = false;

    bool prepared = false;
};

}
