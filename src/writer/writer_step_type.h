#pragma once

namespace fem {
namespace reader {

enum class WriterStepType {
    Static,
    Dynamic,
    Eigenfrequency,
    Buckling
};

} // namespace reader
} // namespace fem
