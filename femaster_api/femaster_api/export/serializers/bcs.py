"""FEMaster serialization for boundary conditions."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword, target_token, trim_missing
from femaster_api.model.supports import SupportCollectorRepository


def write_boundaries(support_collectors: SupportCollectorRepository) -> str:
    blocks: list[str] = []
    for collector in support_collectors.all():
        for support in collector.supports:
            values = trim_missing(support.values)
            blocks.append(
                block(
                    [
                        keyword(
                            "SUPPORT",
                            SUPPORT_COLLECTOR=collector.name,
                            ORIENTATION=None if support.orientation is None else support.orientation.name,
                        ),
                        csv((target_token(support.target), *values)),
                    ]
                )
            )
    return "\n\n".join(block for block in blocks if block)
