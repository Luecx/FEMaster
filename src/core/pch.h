/**
 * @file pch.h
 * @brief Precompiled header for frequently used Eigen and standard-library headers.
 *
 * Keep backend-specific solver headers out of this file. In particular,
 * PardisoSupport and CUDA headers belong only to translation units that use
 * those APIs directly. types_eig.h intentionally remains here so MKL-enabled
 * builds define EIGEN_USE_MKL_ALL before Eigen is parsed by the PCH.
 */

#pragma once

#include "types_eig.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
