import copy
import numpy as np

class LatticeConf:
    """
    Class that stores a lattice system configuration particle coordinates with methods for
    manipulation of the configuration data 

    Note: will need to be refined - divide different trajectory data to
    attributes applicable to the trajectory as a whole and data. This is
    similar to hdf5 format.
    """

    def __init__(self, config, inplace=True):
        """Creates a configuration object from a dictionary or another config
        object"""

        if traj is None:
            self._config = {}

        elif isinstance(traj, dict):
            if inplace:
                self._config = config
            else:
                self._config = copy.deepcopy(config)

        elif isinstance(config, LatticeConf):
            if inplace:
                self._config = config._config
            else:
                self._config = copy.deepcopy(config._config)


    def __getitem__(self, key):
        """
        Implement column (property) access and row slicing of trajectories

        Parameters
        ----------
        key : str or unicode
            if str, select appropriate property (e.g., box)

        Returns
        -------
        data: configuration attribute
        """

        elif isinstance(key, str) or isinstance(key, unicode):
            return self._config[key]

        else:
            raise TypeError('Invalid argument type: {}: {}.'.format(key, type(key)))



    def to_xyz(self, file_name):
        """
        Save configuration to extended XYZ file format
        """

        with open(file_name, 'w') as f:
            # cycle through configurations in trajectory, assign atom names, and write to file

            for box, xyz in zip(self['box'], self['xyz']):

                # write total number of atoms
                nat = sum(self['atom_num'])
                #f.write(f'{nat}\n')
                f.write('{}\n'.format(nat))

                # write box parameters
                ax, ay, az = box[0,0], box[0,1], box[0,2]
                bx, by, bz = box[1,0], box[1,1], box[1,2]
                cx, cy, cz = box[2,0], box[2,1], box[2,2]
                #f.write(f'{ax} {ay} {az} {bx} {by} {bz} {cx} {cy} {cz}\n')
                f.write('{} {} {} {} {} {} {} {} {}\n'.format(ax, ay, az, bx, by, bz, cx, cy, cz))

                # write atom coordinates
                i = 0
                for atom_name, atom_num in zip(self['atom_name'], self['atom_num']):
                    for _ in range(atom_num):
                        x, y, z = (box.T).dot(xyz[i])
                        #f.write(f'{atom_name} {x} {y} {z}\n')
                        f.write('{} {} {} {}\n'.format(atom_name, x, y, z))
                        i += 1

