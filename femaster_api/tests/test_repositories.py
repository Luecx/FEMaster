from femaster_api import Part, PartRepository


def test_named_part_access_and_default_invariant():
    parts = PartRepository()
    wing = parts.add(Part("WING"))
    bolt = parts.add(Part("BOLT"))

    assert parts.default() is parts[0]
    assert parts["WING"] is wing
    assert parts[1] is wing
    assert parts["BOLT"] is bolt
    assert parts[2] is bolt


def test_default_part_cannot_be_removed():
    parts = PartRepository()

    for key in (0, PartRepository.DEFAULT_NAME):
        try:
            parts.remove(key)
        except ValueError:
            pass
        else:
            raise AssertionError("default part removal must fail")


def test_clear_preserves_default_part():
    parts = PartRepository()
    parts.add(Part("WING"))
    parts.clear()

    assert len(parts) == 1
    assert parts.default() is parts[0]
