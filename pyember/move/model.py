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

class MMCMove:
    """Class managing mmc moves"""

    self._supported_moves = set([
                'spin_flip_3d'
            ])

    def __init__(self, move_type, config):
        self.move_type = move_type

        if move_type == 'spin_flip_3d':
            self.move = self.spin_flip_3d_propose
            self.accept = self.spin_flip_3d_accept

        self.boxvec = np.diag(config['box'])


    def spin_flip_3d_propose(self, config):
        """Select a random spin from a given configuration and generate its random orientation"""

        ix = np.random.randint(self.boxvec[0])
        iy = np.random.randint(self.boxvec[1])
        iz = np.random.randint(self.boxvec[2])

        sz = 2*np.random.random() - 1
        st = np.sqrt(1 - sz*sz)
        phi = 2*np.pi*np.random.random()
        sx = st*np.sin(phi)
        sy = st*np.cos(phi)

        return np.array([ix, iy, iz]), np.array([sx, sy, sz])


    def spin_flip_3d_accept(self, config, ri, si, hamilton):
        config['latt_intra'][ri,:] = si

        i = config['latt_i'][ri]
        hamilton.energy_i[i] += hamilton.dui[0] + 0.5*np.sum(hamilton.dui[1:])

        for nbr, du in zip(self.nbrlist, hamilton.dui[1:]):
            rj = tuple((np.array(ri) + nbr) % self.boxvec)
            j = config['latt_i'][rj]
            hamilton.energy_i[j] += 0.5*hamilton.dui[ni]

        hamilton.energy_total += np.sum(hamilton.dui)

