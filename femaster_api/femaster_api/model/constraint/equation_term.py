"""One linear-equation term."""


class EquationTerm:
    """Coefficient multiplying one nodal degree of freedom."""

    def __init__(self, node: int | str, dof: int, coefficient: float) -> None:
        self.node = node
        self.dof = int(dof)
        self.coefficient = float(coefficient)
