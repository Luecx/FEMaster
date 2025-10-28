from __future__ import annotations

from pathlib import Path

if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    __package__ = "fembridge"

from .core import Model, Runner
from .geometry import StraightSegment, SegmentGroup
from .materials import Material, ElasticityIsotropic, Density
from .sections.solid_section import SolidSection
from .sets.elementset import ElementSet
from .supports import Support, SupportCollector
from .loads import CLoad, LoadCollector
from .steps import LinearStaticStep, LinearBucklingStep


def _find_nodeset(model: Model, name: str):
    for node_set in model.node_sets._items:
        if node_set is not None and node_set.name == name:
            return node_set
    raise KeyError(f"NodeSet '{name}' not found in model.")


def build_solid_model() -> tuple[Model, SupportCollector, LoadCollector]:
    seg_bottom = StraightSegment((0.0, 0.0), (10.0, 0.0), n_subdivisions=10, name="BOTTOM")
    seg_right  = StraightSegment((10.0, 0.0), (10.0, 10.0), n_subdivisions=10, name="RIGHT")
    seg_top    = StraightSegment((10.0, 10.0), (0.0, 10.0), n_subdivisions=10, name="TOP")
    seg_left   = StraightSegment((0.0, 10.0), (0.0, 0.0), n_subdivisions=10, name="LEFT")

    outer = SegmentGroup([seg_bottom, seg_right, seg_top, seg_left], name="OUTER")

    plate_model = Model.mesh_2d([outer], mesh_type=2, name="Plate10x10")
    solid_model = plate_model.extruded(n=3, spacing=0.3)
    # solid_model = solid_model.subdivided(1)
    solid_model.name = "Solid10x10x3"

    all_elements = [elem for elem in solid_model.elements._items if elem is not None]
    solid_elset = ElementSet("SOLID", all_elements)
    solid_model.add_elementset(solid_elset)

    steel = solid_model.add_material(
        Material("STEEL")
        .set_elasticity(ElasticityIsotropic(210e9, 0.30))
        .set_density(Density(7850.0))
    )
    solid_model.add_section(SolidSection(material=steel, elset=solid_elset))

    left_nodes = _find_nodeset(solid_model, "LEFT")
    right_nodes = _find_nodeset(solid_model, "RIGHT")

    support = Support(left_nodes, (0, 0, 0, 0, 0, 0))
    support_collector = SupportCollector("SUPPORTS_LEFT")
    support_collector.add(support)
    solid_model.add_supportcollector(support_collector)

    load = CLoad(right_nodes, (-1.0e4, 0.0, 0.0))
    load_collector = LoadCollector("LOAD_RIGHT")
    load_collector.add(load)
    solid_model.add_loadcollector(load_collector)

    return solid_model, support_collector, load_collector


def main() -> None:
    model, support_collector, load_collector = build_solid_model()
    model.add_step(
        LinearStaticStep(
            load_collectors=[load_collector],
            support_collectors=[support_collector],
        )
    )
    model.add_step(
        LinearBucklingStep(
            num_eigenvalues=5,
            load_collectors=[load_collector],
            support_collectors=[support_collector],
        )
    )

    print(model.to_femaster())
    print()
    print(model.steps.to_femaster())

    engine_path = Path(__file__).resolve().parent.parent / "bin" / "FEMaster"
    runner = Runner().set_model(model)
    runner.set_engine(Runner.Engine.FEMASTER, path=engine_path)
    runner.set_option(Runner.Option.NO_TEMP_FILES, True)
    solution = runner.run()
    print()
    print(solution)


if __name__ == "__main__":
    main()
