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

    self._required_keys = set(
                ['latt_i', 'latt_type', 'xyz', 'atom_types',
                 'latt_atoms', 'nat', 'pbc', 'box']
            )

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

        if not self._check_consistency():
            raise ValueError('Configuration is not formed properly')


    def __getitem__(self, key):
        """
        Implements column (property) access and row slicing of trajectories

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


    def __setitem__(self, key, value):
        self._config[key] = value

    def _check_consistency(self):

        config_keys = self._config.keys()

        for key in self._config.keys():
            if key not in self._required_keys:
                return False


        return True

