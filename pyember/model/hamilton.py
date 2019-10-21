#!//anaconda/envs/py36/bin/python
#
# File name:   kmc_pld.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import random
import numpy as np
from itertools import product
from collections import Counter, defaultdict
#import events
from .events import EventTree

class Heisenberg:
    """Class defining the Heisenberg Hamiltonian"""

    def __init__(self, config, J=0.0, H=0.0):
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

        self._check_consistency(config)
        self.latt_type = config['latt_type'].lower()
        self.J = J
        self.H = H
        self._setup_neighbors()
        self.energy_i = np.zeros(config['atom_type'].shape, dtype=np.float64)
        self.energy_total = 0.0
        self.boxvec = np.diag(config['box'])
        assert len(self.boxvec.shape) == 1, "Lattice box dimensions are not a vector"


    def _setup_neighbors(self):
        """Creates lists of neighboring sites"""

        nbrlist = []

        if self.latt_type == 'sc':
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


    def get_energy_diff_i(self, config, ri, skr):
        """Returns interaction energy of atom i
        
        Parameters
        ----------
        config: Config object
        ri: nd.array, shape=(3,)
            position of particle i on lattice
        skr: nd.array, shape(3,)
            trial spin coordinates of particle i

        Returns
        -------
        self.dui: np.array, shape(len(nbr)+1,)
            contribution of external and interspin interactions to energy
        """

        sir = config['latt_intra'][ri,:]

        self.dui[0] = H*sir

        for j, nbr in enumerate(self.nbrlist, 1):
            rj = tuple((np.array(ri) + nbr) % self.boxvec)
            sjr = config['latt_intra'][ri,:]
            self.dui[j] = J*(skr.dot(sjr) - sir.dot(sjr)

        return self.dui

    def get_energy_i(self, config, ri):
        """Returns interaction energy of atom i"""

        sir = config['latt_intra'][ri,:]

        self.ui[0] = H*sir

        for j, nbr in enumerate(self.nbrlist, 1):
            rj = tuple((np.array(ri) + nbr) % self.boxvec)
            sjr = config['latt_intra'][ri,:]
            self.ui[j] = J*sir.dot(sjr)

        return self.ui


    def get_energy_total(self, config):
        """Returns interaction energy of the whole lattice system"""

        u_tot = 0.0
        for ri in config['xyz']:
            ui = self.energy(config, ri)
            u_tot += 2*ui[0]
            u_tot += np.sum(ui[1:])

        self.energy_total = 0.5*u_tot
        return self.energy_total


    def _check_consistency(self, config):
        """Checks if the configuration has two internal spin coordinates"""

        if 'latt_intra' not in config.keys():
            raise KeyError('latt_intra key not in configuration for Heisenberg model')
        

