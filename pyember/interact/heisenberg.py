import numpy as np

class Heisenberg:
    """Class defining the Heisenberg Hamiltonian"""

    def __init__(self, config, params):
        """Initializes Heisenberg Hamiltonian

        Parameters
        ----------
        latt_type: str
            lattice type currently supported only simple cubic, SC
        J: float
            spin-spin interaction parameter
        H: float
            spin-external field interaction parameter
        """

        self.latt_type = 'SC_n3'
        self._check_consistency(config)
        self.temp = params['Temp']
        self.beta = 1.0/self.temp
        self.J = params['J']
        self.H = np.array(params['H'])
        self._setup_neighbors()
        self.energy_i = np.zeros(config['atom_types'].shape, dtype=np.float64)
        self.energy_total = 0.0
        self.boxvec = np.diag(config['box'])
        assert len(self.boxvec.shape) == 1, "Lattice box dimensions are not a vector"


    def _setup_neighbors(self):
        """Creates lists of neighboring sites"""

        nbrlist = []

        if 'sc' in self.latt_type.lower():
            # NN
            nbrlist.append(np.array([ 1, 0, 0]))
            nbrlist.append(np.array([-1, 0, 0]))
            nbrlist.append(np.array([ 0, 1, 0]))
            nbrlist.append(np.array([ 0,-1, 0]))
            nbrlist.append(np.array([ 0, 0, 1]))
            nbrlist.append(np.array([ 0, 0,-1]))
        else:
            raise ValueError(f'Chosen {self.latt_type} lattice. Currently only SC lattice is supported for Heisenberg model.')

        self.nbrlist = nbrlist
        self.ui = np.zeros(len(nbrlist) + 1, dtype=np.float64)
        self.dui = np.zeros(len(nbrlist) + 1, dtype=np.float64)


    def get_energy_diff_i(self, config, event):
        """Returns interaction energy of atom i
        
        Parameters
        ----------
        config: Config object
        event: tuple
            Event: initial to final state

        Returns
        -------
        self.dui: np.array, shape(len(nbr)+1,)
            contribution of external and interspin interactions to energy
        """

        ri = event[0][0]
        ix, iy, iz = ri
        sir = event[0][1]
        skr = event[1][1]

        sir = config['latt_intra'][ix, iy, iz, :]

        self.dui[0] = self.H.dot(sir)

        for j, nbr in enumerate(self.nbrlist, 1):
            jx, jy, jz = (np.array(ri) + nbr) % self.boxvec
            sjr = config['latt_intra'][jx, jy, jz, :]

            self.dui[j]  = skr[0]*sjr[0] + skr[1]*sjr[1] + skr[2]*sjr[2]
            self.dui[j] -= sir[0]*sjr[0] + sir[1]*sjr[1] + sir[2]*sjr[2]
            self.dui[j] *= self.J

        return np.sum(self.beta*self.dui)


    def get_energy_i(self, config, ri):
        """Returns interaction energy of atom i"""

        ix, iy, iz = ri
        sir = config['latt_intra'][ix, iy, iz, :]

        self.ui[0] = self.H.dot(sir)

        for j, nbr in enumerate(self.nbrlist, 1):
            jx, jy, jz = (np.array(ri) + nbr) % self.boxvec
            sjr = config['latt_intra'][jx, jy, jz, :]
            self.ui[j] = self.J*sir.dot(sjr)

        return self.ui


    def get_energy_total(self, config):
        """Returns interaction energy of the whole lattice system"""

        u_tot = 0.0
        for ri in config['xyz']:
            ui = self.get_energy_i(config, ri)
            u_tot += 2*ui[0]
            u_tot += np.sum(ui[1:])

        self.energy_total = 0.5*u_tot
        return self.energy_total


    def _check_consistency(self, config):
        """Checks if the configuration has two internal spin coordinates"""

        assert config['latt_type'] == self.latt_type, f"Incompatible lattice {config['latt_type']} vs {self.latt_type}"

        if 'latt_intra' not in config.keys():
            raise KeyError('latt_intra key not in configuration for Heisenberg model')
        
    def get_magnetization(self, config):

        sum_sx = np.sum(config['latt_intra'][:,:,:,0])
        sum_sy = np.sum(config['latt_intra'][:,:,:,1])
        sum_sz = np.sum(config['latt_intra'][:,:,:,2])

        mag = np.sqrt(sum_sx**2 + sum_sy**2 + sum_sz**2)

        return mag, sum_sx, sum_sy, sum_sz
