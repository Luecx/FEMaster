#pragma once

#include "section.h"

namespace fem {

struct TrussSection : Section {
    using Ptr = std::shared_ptr<TrussSection>;

    Precision A = Precision(0);

    void info() override;
    std::string str() const override;
};

} // namespace fem
