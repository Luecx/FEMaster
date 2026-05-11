from femaster_api import LoadCollector, Model, NodalForce, Node, NodeSet, Support, SupportCollector


def test_support_collector_owns_inline_supports() -> None:
    model = Model("collector_support")
    node = model.nodes.add(Node(0.0, 0.0, 0.0))
    fixed = model.sets.add(NodeSet("FIXED", (node,)))

    supports = SupportCollector("BCS").add(
        Support(
            target=fixed,
            values=(0, 0, 0, 0, 0, 0),
        )
    )
    model.support_collectors.add(supports)

    assert model.validate().errors == ()


def test_load_collector_inline_load_targets_are_validated() -> None:
    model = Model("collector_load")
    node = model.nodes.add(Node(0.0, 0.0, 0.0))
    foreign_set = NodeSet("FOREIGN", (node,))

    loads = LoadCollector("LOADS").add(
        NodalForce(
            target=foreign_set,
            values=(1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        )
    )
    model.load_collectors.add(loads)

    errors = model.validate().errors
    assert len(errors) == 1
    assert errors[0].location == "loads"
    assert errors[0].message == "target references a set outside this model"
