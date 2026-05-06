"""FEMaster serialization for sections and beam profiles."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.sections import BeamSection, SectionRepository, ShellSection, SolidSection, TrussSection


def write_sections(sections: SectionRepository) -> str:
    blocks: list[str] = []
    for profile in sections.profiles():
        blocks.append(
            block(
                [
                    keyword("PROFILE", NAME=profile.name),
                    csv((profile.area, profile.iy, profile.iz, profile.j, profile.iyz, profile.ey, profile.ez, profile.refy, profile.refz)),
                ]
            )
        )

    for section in sections.all():
        if isinstance(section, SolidSection):
            blocks.append(keyword("SOLIDSECTION", ELSET=section.element_set.name, MATERIAL=section.material.name, ORIENTATION=_name(section.orientation)))
        elif isinstance(section, ShellSection):
            blocks.append(
                block(
                    [
                        keyword("SHELLSECTION", ELSET=section.element_set.name, MATERIAL=section.material.name, ORIENTATION=_name(section.orientation)),
                        csv((section.thickness,)),
                    ]
                )
            )
        elif isinstance(section, BeamSection):
            lines = [keyword("BEAMSECTION", ELSET=section.element_set.name, MATERIAL=section.material.name, PROFILE=section.profile.name)]
            if section.orientation is not None:
                lines.append(csv(section.orientation))
            blocks.append(block(lines))
        elif isinstance(section, TrussSection):
            blocks.append(
                block(
                    [
                        keyword("TRUSSSECTION", ELSET=section.element_set.name, MATERIAL=section.material.name),
                        csv((section.area,)),
                    ]
                )
            )
    return "\n\n".join(block for block in blocks if block)


def _name(value):
    return None if value is None else value.name
