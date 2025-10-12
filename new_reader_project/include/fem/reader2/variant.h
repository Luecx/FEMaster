
#pragma once
#include <vector>
#include "segment.h"
#include "condition.h"

namespace fem { namespace reader2 {
struct Variant {
    CondPtr condition_;
    std::vector<Segment> segments_;
    static Variant make(){ return {}; }
    Variant& when(CondPtr c){ condition_ = std::move(c); return *this; }
    Variant& segment(Segment s){ segments_.push_back(std::move(s)); return *this; }
};
}} // ns
