#pragma once

#include "../model/model.h"

namespace fem{
namespace loadcase{

// every load case has loads and supports
struct LoadCase{

    virtual void run() = 0;

};


}
}