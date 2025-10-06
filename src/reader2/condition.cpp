//
// Created by Luecx on 06.10.2025.
//

#include "condition.h"

namespace fem::reader2 {
bool Condition::operator()(const Keyword &kw) const {
    return test ? test(kw) : false;
}
}