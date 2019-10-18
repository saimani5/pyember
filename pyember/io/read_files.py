import numpy as np

def read_xyz(filename):
    """
    Reads extended xyz file with a single configuration

    Parameters
    ----------
    filename: str
              full path and name of the xyz file

    Returns
    -------
    config: dict
          trajectory information with keys given below
    box: list of 3x3 ndarrays
          Box dimensions
    xyz: list of natom x 3 ndarrays
          Atomic configurations
    atom_name: list of str
                atom types (names)
    atom_num: list of ints
                atom numbers for each type
    """

    config = {}

    with open(filename, 'r') as f:
        # number of atoms (spins)
        config['nat'] = int(re.findall('\S+', f.readline())[0])

        # box parameters (type, dimension, shape, periodicity)
        sarr = re.findall('\S+', f.readline())
        dims = tuple(map(int, sarr[1:4]))
        config['box'] = np.diag(dims)
        config['pbc'] = tuple(map(int, sarr[4:7]))
        if len(sarr) > 7:
            dim_intra = len(sarr) - 7

        atom_types = []
        xyz = []
        config['latt_types'] = np.zeros(dims, dtype=int)
        config['latt_intra'] = np.zeros(tuple(dims) + (dim_intra,), dtype='float64')
        for i in range(nat):
            sarr = re.findall('\S+', f.readline())
            t = int(sarr[0])
            r = tuple(map(float, sarr[1:4]))

            atom_types.append(t)
            xyz.append(r)

            config['latt_types'][r] = t
            if dim_intra > 0:
                ir = tuple(map(float, sarr[4:4+dim_intra]))
                config['latt_intra'][r,...] = ir

    config['atom_types'] = np.array(atom_types)
    config['xyz'] = np.array(xyz)
        
    return config
