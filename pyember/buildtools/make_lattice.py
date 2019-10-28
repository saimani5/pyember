#!//anaconda/envs/py36/bin/python
#
# File name:   mklatt.py
# Date:        2018/08/02 16:15
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
from itertools import product
import numpy as np

def make_sc(box):
    print('sc', *box)
    return

def make_bcc(box):
    print('bcc', *box)
    return

def make_fcc(box):

    lx, ly, lz = box

    latt = {}
    latt['nat'] = lx*ly*lz
    latt['box'] = ['fcc', 2*lx, ly, lz]
    latt['xyz'] = []

    # box dimensions in lattice units

    # layer number
    for iz in range(lz):
        # layer structure
        for iy in range(ly):
            for ix in range(lx):
                rx = 2*ix + (iy + iz)%2
                latt['xyz'].append(['Ni', rx, iy, iz])

    return latt

def make_heisenberg(dims=(8, 8, 8), pbc=(1, 1, 1), random=True):

    config = {}
    config['nat'] = np.prod(dims)
    config['box'] = np.diag(dims)
    config['pbc'] = pbc

    config['atom_types'] = np.zeros((config['nat'],), dtype=int)
    xyz = [(x, y, z) for x, y, z in product(range(dims[0]), range(dims[1]), range(dims[2]))]
    config['xyz'] = np.array(xyz)

    config['latt_type'] = 'SC_n3'
    config['latt_atoms'] = np.zeros(dims, dtype=int)
    config['latt_intra'] = np.zeros(tuple(dims) + (3,), dtype='float64')

    if random:
        # new spin orientation
        sz = 2*np.random.random(size=dims) - 1
        st = np.sqrt(1 - sz*sz)
        phi = 2*np.pi*np.random.random(size=dims)
        sx = st*np.sin(phi)
        sy = st*np.cos(phi)

        config['latt_intra'][...,0] = sx
        config['latt_intra'][...,1] = sy
        config['latt_intra'][...,2] = sz

        r2 = sx**2 + sy**2 + sz**2 - 1.0
        assert np.max(r2) < 1e-8, "Spin magnitude not close to 1"

    return config


if __name__ == "__main__":

    # dictionary of lattice build functions
    make_lattice = {'fcc':make_fcc, 'bcc':make_bcc, 'sc':make_sc}

    # read lattice type and parameters
    with open(sys.argv[1], 'r') as f:

        # lattice type
        ltype = re.findall('\S+', f.readline())[0]

        # box dimensions
        box = [int(d) for d in re.findall('\S+', f.readline())]

    # make lattice
    latt = make_lattice[ltype](box)

# end of mklatt.py 
