
#pragma once
#include <vector>
#include <memory>
#include "pattern_element.h"

namespace fem { namespace reader2 {
struct Pattern {
    bool multiline = false;
    std::vector<std::shared_ptr<PatternElementBase>> elems;
    static Pattern make(){ return {}; }
    Pattern& allow_multiline(){ multiline = true; return *this; }
    template<class T, std::size_t N>
    Pattern& fixed(){ elems.emplace_back(std::make_shared<Fixed<T,N>>()); return *this; }
};
}} // ns
