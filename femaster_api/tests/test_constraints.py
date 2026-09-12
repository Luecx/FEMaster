"""Constraint API checks for nested helper types and object-valued relations.

Constraint helpers whose lifecycle belongs entirely to one constraint remain
nested (`Equation.Term`, `Connector.Type`).  Cross-model relations are different:
those always hold the actual ``Node``, ``NodeRegion`` or ``CoordinateSystem``
objects and explicitly reject textual substitutes.
"""

import pytest

import femaster_api
from femaster_api import (
    Connector,
    Equation,
    Node,
    NodeRegion,
    RectangularCoordinateSystem,
)


def test_equation_term_is_nested_and_holds_node_object():
    node_a = Node(17, 0.0, 0.0, 0.0)
    node_b = Node(42, 1.0, 0.0, 0.0)
    term = Equation.Term(node_a, 2, -1.5)
    equation = Equation((term, (node_b, 1, 3.0)))

    assert len(equation.terms) == 2
    assert isinstance(equation.terms[0], Equation.Term)
    assert equation.terms[0].node is node_a
    assert equation.terms[0].dof == 2
    assert equation.terms[0].coefficient == -1.5

    text = equation.export()
    assert "*EQUATION" in text
    assert "17" in text
    assert "42" in text
    assert not hasattr(femaster_api, "EquationTerm")

    with pytest.raises(TypeError):
        Equation.Term("bolt.17", 2, -1.5)


def test_connector_type_is_nested_and_relations_are_objects():
    master_node = Node(1, 0.0, 0.0, 0.0)
    slave_node = Node(2, 1.0, 0.0, 0.0)
    master = NodeRegion("MASTER", (master_node,))
    slave = NodeRegion("SLAVE", (slave_node,))
    orientation = RectangularCoordinateSystem(
        "GLOBAL",
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
    )

    connector = Connector(
        Connector.Type.RIGID,
        master,
        slave,
        orientation,
    )

    assert connector.type is Connector.Type.RIGID
    assert connector.nset1 is master
    assert connector.nset2 is slave
    assert connector.coordinate_system is orientation
    assert "TYPE=RIGID" in connector.export()
    assert "NSET1=MASTER" in connector.export()
    assert not hasattr(femaster_api, "ConnectorType")

    with pytest.raises(TypeError):
        Connector(Connector.Type.RIGID, "MASTER", slave, orientation)

    with pytest.raises(TypeError):
        Connector("RIGID", master, slave, orientation)
