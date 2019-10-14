"""The configurations class should have

1. 3D Array
2. List of atoms and coordinates
3. Conversion between 1 and 2
4. Periodic boundary conditions
5. Lattice Type
6. Atom Names
7. Helper functions (possibly)
8. Easier slicing methods
"""

import numpy as np

class SpinCubic(Universe):
    """Condiguration space of the classical n-vector models"""

    def __init__(self, incfg_file=None, rand_init=True, spatial_dim=(8, 8, 8), boundary_cond=(1, 1, 1), n_vector=1):

        if incfg_file in not None:
            # configuration xyz file should contain all that is necessary
            self.config, self.lattice_dim, self.boundary_cond = read_xyzfile(incfg_file)

        else:

            assert len(lattice_dim) == len(boundary_cond), "Mismatched lattice/boundary dimensions"

            self.spatial_dim = tuple(spatial_dim)
            self.boundary_cond = tuple(boundary_cond)

            assert n_vector > 0 and n_vector < 4, f"Unsupported n_vector dimension: {n_vector}."
            self.n_vector = n_vector

            self.lattice_dim = self.spatial_dim
            if n_vector > 1: 
                self.lattice_dim += (n_vector - 1, )

            # create an simulation box, each site with spin coordinates set to zero
            self.config = np.zeros(self.lattice_dim, dtype=np.float64)

            # if random configuration initialization, set spins of length 1 to random directions
            if rand_init:
                self.config = np.random


