"""FEMaster serialization for loads and amplitudes."""

from __future__ import annotations

from femaster_api.export.context import ExportContext
from femaster_api.export.femaster_format import block, csv, keyword, trim_missing
from femaster_api.model.loads import (
    InertialLoad,
    LoadCollectorRepository,
    LoadRepository,
    NodalForce,
    PressureLoad,
    SurfaceTraction,
    ThermalLoad,
    VolumeLoad,
)


def write_loads(loads: LoadRepository, load_collectors: LoadCollectorRepository, context: ExportContext) -> str:
    blocks: list[str] = []
    for amplitude in loads.amplitudes():
        lines = [keyword("AMPLITUDE", NAME=amplitude.name, TYPE=amplitude.interpolation.value)]
        lines.extend(csv(point) for point in amplitude.points)
        blocks.append(block(lines))

    for collector in load_collectors:
        for entry in collector.loads:
            if isinstance(entry, NodalForce):
                blocks.append(
                    block(
                        [
                            keyword(
                                "CLOAD",
                                LOAD_COLLECTOR=collector.name,
                                ORIENTATION=_name(entry.orientation),
                                AMPLITUDE=_name(entry.amplitude),
                            ),
                            csv((context.target_token(entry.target), *trim_missing(entry.values))),
                        ]
                    )
                )
            elif isinstance(entry, SurfaceTraction):
                blocks.append(
                    block(
                        [
                            keyword(
                                "DLOAD",
                                LOAD_COLLECTOR=collector.name,
                                ORIENTATION=_name(entry.orientation),
                                AMPLITUDE=_name(entry.amplitude),
                            ),
                            csv((context.target_token(entry.target), *entry.values)),
                        ]
                    )
                )
            elif isinstance(entry, PressureLoad):
                blocks.append(
                    block(
                        [
                            keyword("PLOAD", LOAD_COLLECTOR=collector.name, AMPLITUDE=_name(entry.amplitude)),
                            csv((context.target_token(entry.target), entry.pressure)),
                        ]
                    )
                )
            elif isinstance(entry, VolumeLoad):
                blocks.append(
                    block(
                        [
                            keyword(
                                "VLOAD",
                                LOAD_COLLECTOR=collector.name,
                                ORIENTATION=_name(entry.orientation),
                                AMPLITUDE=_name(entry.amplitude),
                            ),
                            csv((context.target_token(entry.target), *entry.values)),
                        ]
                    )
                )
            elif isinstance(entry, ThermalLoad):
                blocks.append(
                    keyword(
                        "TLOAD",
                        LOAD_COLLECTOR=collector.name,
                        TEMPERATUREFIELD=entry.temperature_field.name,
                        REFERENCETEMPERATURE=entry.reference_temperature,
                    )
                )
            elif isinstance(entry, InertialLoad):
                blocks.append(
                    block(
                        [
                            keyword(
                                "INERTIALOAD",
                                LOAD_COLLECTOR=collector.name,
                                CONSIDER_POINT_MASSES=int(entry.consider_point_masses),
                            ),
                            csv(
                                (
                                    entry.target.name,
                                    *entry.center,
                                    *entry.center_acceleration,
                                    *entry.omega,
                                    *entry.alpha,
                                )
                            ),
                        ]
                    )
                )
    return "\n\n".join(block for block in blocks if block)


def _name(value):
    return None if value is None else value.name
