
class Material:
    """Base material class."""
    def __init__(self, name: str, mtype: str):
        if not name:
            raise ValueError("Material name must be non-empty.")
        self.name = str(name)
        self.mtype = str(mtype)

    def to_femaster(self) -> str:
        raise NotImplementedError("Material.to_femaster must be implemented by subclasses.")

    def to_asami(self) -> str:
        raise NotImplementedError("Material.to_asami not implemented yet.")
