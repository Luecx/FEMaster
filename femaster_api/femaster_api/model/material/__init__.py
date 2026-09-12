"""Global FEMaster material and elasticity definitions.

A project owns named ``Material`` objects, while each constitutive elasticity
form is modeled by its own class and module.  Sections reference materials by
semantic name, so constitutive definitions stay globally reusable and remain
independent of part-local element assignment.
"""

from .material import Material
from .material_elasticity import Elasticity
from .material_elasticity_abd import ABDElasticity
from .material_elasticity_generalized_isotropic import GeneralizedIsotropicElasticity
from .material_elasticity_isotropic import IsotropicElasticity
from .material_elasticity_orthotropic import OrthotropicElasticity
from .material_repository import MaterialRepository

__all__ = [
    "ABDElasticity",
    "Elasticity",
    "GeneralizedIsotropicElasticity",
    "IsotropicElasticity",
    "Material",
    "MaterialRepository",
    "OrthotropicElasticity",
]
