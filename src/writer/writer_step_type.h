#pragma once

namespace fem {
namespace io {
namespace writer {

enum class WriterStepType {
    Static,
    Dynamic,
    Eigenfrequency,
    Buckling
};

} // namespace writer
} // namespace io
} // namespace fem
