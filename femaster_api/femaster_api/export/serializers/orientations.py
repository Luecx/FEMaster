"""FEMaster serialization for coordinate systems."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.orientations import CylindricalOrientation, OrientationRepository, RectangularOrientation


def write_orientations(orientations: OrientationRepository) -> str:
    blocks: list[str] = []
    for orientation in orientations.all():
        if isinstance(orientation, RectangularOrientation):
            data = [*orientation.x_axis]
            if orientation.y_axis is not None:
                data.extend(orientation.y_axis)
            if orientation.z_axis is not None:
                data.extend(orientation.z_axis)
            blocks.append(block([keyword("ORIENTATION", NAME=orientation.name, TYPE="RECTANGULAR"), csv(data)]))
        elif isinstance(orientation, CylindricalOrientation):
            blocks.append(
                block(
                    [
                        keyword("ORIENTATION", NAME=orientation.name, TYPE="CYLINDRICAL"),
                        csv((*orientation.origin, *orientation.axis, *orientation.reference)),
                    ]
                )
            )
    return "\n\n".join(block for block in blocks if block)
