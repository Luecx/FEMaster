//
// Simple printable interface for one-line string representations
//

#pragma once

#include <string>

namespace fem {

struct Printable {
    virtual ~Printable() = default;
    virtual std::string str() const = 0;
};

} // namespace fem

