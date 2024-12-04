//
// Created by f_eggers on 04.12.2024.
//

#ifndef NAMABLE_H
#define NAMABLE_H

#include <string>
#include "../core/core.h"

struct Namable {
    const std::string name;

    Namable(std::string p_name) : name(p_name) {}
};

#endif //NAMABLE_H
