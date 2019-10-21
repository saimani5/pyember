import re
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
        config['latt_type'] = sarr[0]
        dims = tuple(map(int, sarr[1:4]))
        config['box'] = np.diag(dims)
        config['pbc'] = tuple(map(int, sarr[4:7]))
        if len(sarr) > 7:
            dim_intra = len(sarr) - 7

        atom_types = []
        xyz = []
        config['latt_i'] = np.zeros(dims, dtype=int)
        config['latt_atoms'] = np.zeros(dims, dtype=int)
        config['latt_intra'] = np.zeros(tuple(dims) + (dim_intra,), dtype='float64')
        for i in range(config['nat']):
            sarr = re.findall('\S+', f.readline())
            t = int(sarr[0])
            r = tuple(map(int, sarr[1:4]))

            atom_types.append(t)
            xyz.append(r)

            config['latt_i'][r] = i
            config['latt_atoms'][r] = t

            for j in range(dim_intra):
                ci = float(sarr[4 + j])
                config['latt_intra'][r[0], r[1], r[2], j] = ci

    config['atom_types'] = np.array(atom_types)
    config['xyz'] = np.array(xyz)
        
    return config

def write_xyz(config, filename):
    """
    Reads extended xyz file with a single configuration

    Parameters
    ----------
    config: dict or config object to be saved
    filename: str
              full path and name of the xyz file

    """

    with open(filename, 'w') as f:
        # number of atoms (spins)
        f.write("{}\n".format(config['nat']))

        # information line
        f.write("{} ".format(config['latt_type']))
        f.write("{} {} {}".format(*list(np.diag(config['box']))))
        f.write(" {} {} {}".format(*list(config['pbc'])))

        dims_intra = config['latt_intra'].shape[-1]
        for i in range(dims_intra):
            f.write(" 1")

        f.write("\n")

        # coordinates
        for i in range(config['nat']):
            f.write("{} ".format(config['atom_types'][i]))

            ix, iy, iz = list(map(lambda x: int(round(x)), config['xyz'][i]))

            f.write("{} {} {}".format(ix, iy, iz))

            for j in range(dims_intra):
                f.write(" {}".format(config['latt_intra'][ix, iy, iz, j]))

            f.write("\n")

