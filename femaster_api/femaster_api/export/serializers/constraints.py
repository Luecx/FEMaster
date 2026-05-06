"""FEMaster serialization for constraints."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.constraints import (
    ConnectorConstraint,
    ConstraintRepository,
    CouplingConstraint,
    RBMConstraint,
    TieConstraint,
)


def write_constraints(constraints: ConstraintRepository) -> str:
    blocks: list[str] = []
    for item in constraints.all():
        if isinstance(item, RBMConstraint):
            blocks.append(keyword("RBM", ELSET=item.element_set.name))
        elif isinstance(item, CouplingConstraint):
            blocks.append(
                block(
                    [
                        keyword(
                            "COUPLING",
                            MASTER=item.master.name,
                            TYPE=item.type.value,
                            SFSET=item.slave.name if item.slave.entity_type.value == "surface" else None,
                            SLAVE=None if item.slave.entity_type.value == "surface" else item.slave.name,
                        ),
                        csv(tuple(int(value) for value in item.dofs)),
                    ]
                )
            )
        elif isinstance(item, ConnectorConstraint):
            blocks.append(
                keyword(
                    "CONNECTOR",
                    TYPE=item.type.value,
                    NSET1=item.nset1.name,
                    NSET2=item.nset2.name,
                    COORDINATESYSTEM=item.coordinate_system.name,
                )
            )
        elif isinstance(item, TieConstraint):
            blocks.append(
                keyword(
                    "TIE",
                    MASTER=item.master.name,
                    SLAVE=item.slave.name,
                    ADJUST="YES" if item.adjust else "NO",
                    DISTANCE=item.distance,
                )
            )
    return "\n\n".join(block for block in blocks if block)
