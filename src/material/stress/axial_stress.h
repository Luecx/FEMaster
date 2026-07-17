#pragma once

#include "../../core/types_num.h"

namespace fem {

struct AxialStress {
    AxialStress() = default;
    explicit AxialStress(Precision value);

    [[nodiscard]] Precision  value() const;
    [[nodiscard]] Precision& value();

protected:
    Precision value_{};
};

} // namespace fem
