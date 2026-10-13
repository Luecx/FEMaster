#pragma once

#include "profile.h"

#include <string>
#include <vector>

namespace fem {

// Convert supported Abaqus library-section dimensions into FEMaster's
// linear generalized beam properties. Unsupported section geometries must
// be rejected rather than silently using an unrelated approximation.
Profile::Ptr make_beam_geometry_profile(const std::string& name,
                                        const std::string& shape,
                                        const std::vector<Precision>& dimensions);

} // namespace fem
