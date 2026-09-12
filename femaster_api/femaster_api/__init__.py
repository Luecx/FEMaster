"""FEMaster Python API.

The public API mirrors FEMaster's semantic model: Project owns named global
repositories, Part owns local topology, collectors own loads/supports, and each
exportable class writes its own FEMaster representation.
"""

from .fields import Field, FieldDomain, FieldRepository, FieldType
from .io import FrdReader, InpReader, ResReader, read_input, read_result
from .model import *
from .project import Project
from .repository import IdRepository, NamedObject, NamedRepository
from .results import Frame, LoadCase, Result

__all__ = [name for name in globals() if not name.startswith("_")]
