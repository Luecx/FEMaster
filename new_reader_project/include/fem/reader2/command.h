
#pragma once
#include <string>
#include <vector>
#include "variant.h"

namespace fem { namespace reader2 {
struct Command {
    std::string name_;
    CondPtr admit_;
    std::vector<Variant> variants_;
    explicit Command(std::string n): name_(std::move(n)) {}
    Command& allow_if(CondPtr c){ admit_ = std::move(c); return *this; }
    Command& variant(Variant v){ variants_.push_back(std::move(v)); return *this; }
};
}} // ns
