import sys
import os
import random
import copy
import yaml
from .mmcsim import MMCSim
from .kmcsim import KMCSim

class Simulation:
    """ Class for setting up a simulation based on given control information.
    
    Control information is given either in config.yml configuration file or a dict.
    """


    def __init__(self, setup_info=None):
        """
        Initializes simulation object either from an input file or a
        dictionary with appropriate parameters.

        Parameters
        ----------
        setup_info: str or dict
            contains information for setting up a simulation model and flow control
            variables, which are then stored in a sim_parms dict.
        """

        # Read and check basic consistency of the control file/dict
        if isinstance(setup_info, dict):
            self.sim_params = self._check_control_dict(setup_info)
        elif isinstance(setup_info, str):
            with open(setup_info, 'r') as f: 
                setup_dict = yaml.safe_load(f)
            self.sim_params = self._check_control_dict(setup_dict)
        else:
            raise TypeError("Simulation setup info must be either a file name or dict")

        # Build simulation object while checking compatibility of different objects
        if self.sim_params['sim_type'] == 'MMC':
            mmc = MMCSim(self.sim_params)
            self.run = mmc.run
        elif self.sim_params['sim_type'] == 'KMC':
            kmc = KMCSim(self.sim_params)
            self.run = kmc.run


    def _check_control_dict(self, setup_dict):
        """
        Verifies that the supplied dictionary contains the necessary parameters
        and no unknown parameters. Otherwise raises an error.

        Parameters
        ----------
        setup_dict: dict
            Supplied dict of parameter values
        """

        required_pars = set('sim_type', 'config', 'moves', 'time_control')

        # check if all necessary parameters are present
        for par in required_params:
            if par not in setup_dict:
                raise KeyError(f"{par} is required, but not present in the supplied parameters")

        # check if KMC has appropriate rate/barrier information
        if setup_dict['sim_type'] == 'KMC':
            if 'Hamilton' in setup_dict:
                assert 'barriers' in setup_dict, "Barriers are missing for KMC"
            else:
                assert 'rates' in setup_dict, "No rates information for KMC"

        # check basic time control parameters
        assert 'total' in setup_dict['time_control'], "Total simulation length is missing"

        # if missing set non-essential time control parameters to defaults
        if 'save' not in setup_dict['time_control']:
            setup_dict['time_control']['save'] = 10

        if 'print' not in setup_dict['time_control']:
            setup_dict['time_control']['print'] = 10

        if 'measure' not in setup_dict['time_control']:
            setup_dict['time_control']['measure'] = 10


if __name__ == "__main__":

    # start simulation with parameters defined in an control file
    sim = Simulation(sys.argv[0])

    sim.run()

