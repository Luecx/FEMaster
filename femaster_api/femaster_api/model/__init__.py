"""Public FEMaster model and result object graph.

All domain classes live below ``model``.  The package is intentionally organized
by FEM concept—nodes, elements, surfaces, regions, materials, sections, loads,
steps, fields and results—rather than by serializer/importer infrastructure.
``Project`` is the editable model root; ``Result`` is the post-processing root.
"""

from .amplitude import *
from .common import *
from .constraint import *
from .coordinate_system import *
from .element import *
from .feature import *
from .field import *
from .instance import *
from .load import *
from .material import *
from .node import *
from .part import *
from .profile import *
from .project import Project
from .region import *
from .result import *
from .section import *
from .step import *
from .support import *
from .surface import *

__all__ = [name for name in globals() if not name.startswith("_")]
