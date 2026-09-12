"""Constraint API checks for nested helper types and native export semantics.

The public model should expose only actual FEMaster constraint concepts at package
scope.  Helper types whose meaning is entirely owned by one constraint remain
nested on that constraint class.  These tests lock that API shape down while also
checking that the nested representations still export the expected native tokens.
"""

import femaster_api
from femaster_api import Connector, Equation


def test_equation_term_is_nested_and_preserves_named_values():
    term = Equation.Term("bolt.17", 2, -1.5)
    equation = Equation((term, (42, 1, 3.0)))

    assert len(equation.terms) == 2
    assert isinstance(equation.terms[0], Equation.Term)
    assert equation.terms[0].node == "bolt.17"
    assert equation.terms[0].dof == 2
    assert equation.terms[0].coefficient == -1.5

    text = equation.export()
    assert "*EQUATION" in text
    assert "bolt.17" in text
    assert not hasattr(femaster_api, "EquationTerm")


def test_connector_type_is_nested_and_known_strings_are_normalized():
    connector = Connector(
        Connector.Type.RIGID,
        "MASTER",
        "SLAVE",
        "GLOBAL",
    )
    parsed_like_connector = Connector(
        "CARTESIAN",
        "MASTER",
        "SLAVE",
        "GLOBAL",
    )

    assert connector.type is Connector.Type.RIGID
    assert parsed_like_connector.type is Connector.Type.CARTESIAN
    assert "TYPE=RIGID" in connector.export()
    assert not hasattr(femaster_api, "ConnectorType")
