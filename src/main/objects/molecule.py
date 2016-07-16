class Molecule:
    """Molecule object containing nuclei and the molecules point group.

    Attributes
    ----------
    nuclei_array : List[Nuclei]
    point_group : PointGroup

    """
    def __init__(self, nuclei_array, point_group):
        self.nuclei_array = nuclei_array
        self.point_group = point_group
