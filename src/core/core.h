#pragma once

#include "assert.h"
#include "logging.h"
#include "types_eig.h"
#include "timer.h"
#include "config.h"

#ifdef SUPPORT_GPU
#include "../cuda/assert_cuda.h"
#include "../cuda/cuda.h"
#endif