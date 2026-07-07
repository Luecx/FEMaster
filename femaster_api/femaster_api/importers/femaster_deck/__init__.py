"""FEMaster input deck importer."""

from .parser import FEMasterInputError, load_model_from_inp

__all__ = ["FEMasterInputError", "load_model_from_inp"]
