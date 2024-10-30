#pragma once

#include <iostream>
#include <memory>

namespace fem::model {

using ElementType = int;
using ElementTypeFlags = int;

enum ElementTypes : ElementType {
    StructuralType = 1 << 0, // 0001
    ThermalType    = 1 << 1, // 0010
    FluidType      = 1 << 2, // 0100
    N_Types        = 3
};

}

