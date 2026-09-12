"""Assembly-level FEMaster constraint definitions.

Each concrete constraint owns its FEMaster representation. ConstraintRepository
only preserves model order and delegates export; it does not inspect concrete
constraint types or duplicate serializer logic.
"""

from __future__ import annotations

from enum import Enum
from typing import Iterable, Iterator

from .._format import block, csv, keyword


class Constraint:
    """Base class for one assembly-level constraint."""

    def to_femaster(self) -> str:
        raise NotImplementedError


class CouplingType(Enum):
    """Supported FEMaster coupling formulations."""

    KINEMATIC    = "KINEMATIC"
    DISTRIBUTING = "DISTRIBUTING"


class ConnectorType(Enum):
    """Connector kinematic type token passed directly to FEMaster."""

    RIGID = "RIGID"
    CARTESIAN = "CARTESIAN"


class RigidBodyConstraint(Constraint):
    """Rigid-body removal/constraint acting on an element region."""

    def __init__(self, element_region: str) -> None:
        self.element_region = str(element_region)

    def to_femaster(self) -> str:
        return keyword("RBM", ELSET=self.element_region)


class Coupling(Constraint):
    """Kinematic or distributing coupling between master and slave regions."""

    def __init__(
        self,
        master: str,
        slave: str,
        *,
        type: CouplingType = CouplingType.KINEMATIC,
        dofs: Iterable[bool | int] = (1, 1, 1, 1, 1, 1),
        slave_is_surface: bool = False,
    ) -> None:
        self.master           = str(master)
        self.slave            = str(slave)
        self.type             = type
        self.dofs             = tuple(int(bool(value)) for value in dofs)
        self.slave_is_surface = bool(slave_is_surface)

    def to_femaster(self) -> str:
        return block([
            keyword(
                "COUPLING",
                MASTER=self.master,
                TYPE=self.type.value,
                SFSET=self.slave if self.slave_is_surface else None,
                SLAVE=None if self.slave_is_surface else self.slave,
            ),
            csv(self.dofs),
        ])


class Connector(Constraint):
    """Connector relation between two node regions in a local coordinate system."""

    def __init__(self, type: ConnectorType | str, nset1: str, nset2: str, coordinate_system: str) -> None:
        self.type              = type.value if isinstance(type, ConnectorType) else str(type)
        self.nset1             = str(nset1)
        self.nset2             = str(nset2)
        self.coordinate_system = str(coordinate_system)

    def to_femaster(self) -> str:
        return keyword(
            "CONNECTOR",
            TYPE=self.type,
            NSET1=self.nset1,
            NSET2=self.nset2,
            COORDINATESYSTEM=self.coordinate_system,
        )


class Tie(Constraint):
    """Tie a slave node/surface region to a master surface or line region."""

    def __init__(
        self,
        master: str,
        slave: str,
        *,
        adjust: bool = True,
        distance: float | None = None,
    ) -> None:
        self.master   = str(master)
        self.slave    = str(slave)
        self.adjust   = bool(adjust)
        self.distance = None if distance is None else float(distance)

    def to_femaster(self) -> str:
        return keyword(
            "TIE",
            MASTER=self.master,
            SLAVE=self.slave,
            ADJUST="YES" if self.adjust else "NO",
            DISTANCE=self.distance,
        )


class EquationTerm:
    """One coefficient multiplying one nodal degree of freedom."""

    def __init__(self, node: int | str, dof: int, coefficient: float) -> None:
        self.node        = node
        self.dof         = int(dof)
        self.coefficient = float(coefficient)


class Equation(Constraint):
    """Linear multi-point equation containing an arbitrary number of terms."""

    def __init__(self, terms: Iterable[EquationTerm] = ()) -> None:
        self.terms = list(terms)

    def add(self, node: int | str, dof: int, coefficient: float) -> "Equation":
        """Append one equation term and return this equation."""

        self.terms.append(EquationTerm(node, dof, coefficient))
        return self

    def to_femaster(self) -> str:
        values: list[object] = []
        for term in self.terms:
            values.extend((term.node, term.dof, term.coefficient))
        return block([
            keyword("EQUATION"),
            csv((len(self.terms),)),
            csv(values),
        ])


class ConstraintRepository:
    """Ordered collection of heterogeneous assembly-level constraints."""

    def __init__(self) -> None:
        self._items: list[Constraint] = []

    def add(self, constraint: Constraint) -> Constraint:
        if not isinstance(constraint, Constraint):
            raise TypeError("constraint must derive from Constraint")
        self._items.append(constraint)
        return constraint

    def to_femaster(self) -> str:
        return "\n\n".join(item.to_femaster() for item in self._items)

    def __getitem__(self, index: int) -> Constraint:
        return self._items[index]

    def __iter__(self) -> Iterator[Constraint]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
