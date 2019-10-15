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

class SpinLattice(Configuration):
    """Condiguration space of the classical n-vector models"""

    def __init__(self, incfg_file=None, rand_init=True, spatial_dim=(8, 8, 8), boundary_cond=(1, 1, 1), n_vector=1):

        if incfg_file in not None:
            # configuration xyz file should contain all that is necessary
            self.config, self.lattice_dim, self.boundary_cond, self.n_vector = read_spin_xyzfile(incfg_file)

        else:

            assert len(lattice_dim) == len(boundary_cond), "Mismatched lattice/boundary dimensions"

            self.spatial_dim = tuple(spatial_dim)
            self.boundary_cond = tuple(boundary_cond)

            assert n_vector > 0 and n_vector < 4, f"Unsupported n_vector dimension: {n_vector}."
            self.n_vector = n_vector

            self.lattice_dim = self.spatial_dim
            if n_vector == 3: 
                self.lattice_dim += (n_vector - 1, )

            if rand_init:
                # if random configuration initialization, set spins of length 1 to random directions
                self.config = self._random_conf()
            else:
                # else create an simulation box, each site with spin coordinates set to zero
                self.config = self._uniform_conf()


    def _uniform_conf(self):
        """Generage a configuration with all same spins.

        Ising model set to 1, XY and Heisenberg models set to 0.0.
        """

        if self.n_vector == 1:
            # Simple ising model with all spins 1
            self.config = np.ones(self.lattice_dim, dtype=np.float64)

        elif self.n_vector == 2:
            # XY model with angle theta set to 0.0
            self.config = np.zeros(self.lattice_dim, dtype=np.float64)

        elif self.n_vector == 3:
            # Heisenberg model all theta and phi set to 0.0
            self.config = np.zeros(self.lattice_dim, dtype=np.float64)

        else:
            raise ValueError


    def _random_conf(self):
        """Generage random spin configuration"""

        if self.n_vector == 1:
            # Simple ising model (1, -1)
            self.config = 2.*np.random.randint(2, size=self.config.shape).astype('float64') - 1.

        elif self.n_vector == 2:
            # XY model with angle theta uniformly distributed between -pi and +pi
            self.config = np.pi*(2.*np.random.random(size=self.config.shape) - 1.)

        elif self.n_vector == 3:
            # Heisenberg model
            # cos(theta) uniformly distributed in (-1, 1)
            self.config[:,:,:,0] = 2.*np.random.random(size=self.spatial_dim) - 1.
            # phi uniformly distributed in (0, 2*pi)
            self.config[:,:,:,1] = 2.*np.pi*np.random.random(size=self.spatial_dim) 

        else:
            raise ValueError


    def read_spin_xyzfile(self, filename):

        with open(filename, 'r') as f:
            # number of atoms (spins)
            nat = int(re.findall('\S+', f.readline())[0])

            # box parameters (type, dimension, shape, periodicity)
            sarr = re.findall('\S+', f.readline())
            n_vector = int(sarr[0])
            n_dim = int(sarr[1])
            lattice_dim = tuple(map(int, sarr[2: 2+n_dim]))
            boundary_cond = tuple(map(int, sarr[2+n_dim: 2+2*n_dim]))
            
            config = np.empty
            for i in range(nat):
                sarr = re.findall('\S+', f.readline())
                if n_vector == 1:


        return config, lattice_dim, boundary_cond, n_vector

