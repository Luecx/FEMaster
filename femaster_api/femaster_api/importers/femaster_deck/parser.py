"""FEMaster input deck parser."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, TextIO

from femaster_api.model import Model

from .utils import FEMasterInputError, Header, LineStream, parse_header, read_lines


CommandHandler = Callable[["DeckParser", Header], None]


@dataclass
class DeckParser:
    model: Model
    stream: LineStream
    node_by_deck_id: dict[int, object] = field(default_factory=dict)
    element_by_deck_id: dict[int, object] = field(default_factory=dict)
    surface_by_deck_id: dict[int, object] = field(default_factory=dict)
    current_material: object | None = None

    def __post_init__(self) -> None:
        from .commands_analysis import ANALYSIS_COMMANDS
        from .commands_geometry import GEOMETRY_COMMANDS
        from .commands_loads import LOAD_COMMANDS
        from .commands_properties import PROPERTY_COMMANDS

        self.commands: dict[str, CommandHandler] = {}
        self.commands.update(GEOMETRY_COMMANDS)
        self.commands.update(PROPERTY_COMMANDS)
        self.commands.update(LOAD_COMMANDS)
        self.commands.update(ANALYSIS_COMMANDS)

    def run(self) -> Model:
        while True:
            line = self.stream.peek()
            if line is None:
                return self.model
            if not line.startswith("*"):
                raise FEMasterInputError(f"unexpected data outside command: {line!r}")
            header = parse_header(self.stream.pop() or "")
            if header.keyword.startswith("END"):
                continue
            try:
                handler = self.commands[header.keyword]
            except KeyError as exc:
                raise FEMasterInputError(f"unsupported command '*{header.keyword}'") from exc
            handler(self, header)

    def consume_data_lines(self):
        while True:
            line = self.stream.peek()
            if line is None or line.startswith("*"):
                break
            yield self.stream.pop() or ""

    def next_header_is_top_level(self, loadcase_subcommands: set[str]) -> bool:
        line = self.stream.peek()
        if line is None or not line.startswith("*"):
            return False
        keyword = parse_header(line).keyword
        return keyword not in loadcase_subcommands and keyword in self.commands


def load_model_from_inp(source: str | Path | TextIO | Iterable[str], *, model: Model | None = None) -> Model:
    lines = read_lines(source)
    parser = DeckParser(model or Model(), LineStream(lines))
    return parser.run()
