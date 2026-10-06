#pragma once

namespace fem::requests {

enum NodalRequests {
    DISPLACEMENT,
    REACTION_FORCE,
    STRESS,
    MISES,
    SECTION_FORCE,
    SECTION_MOMENT,
};

}