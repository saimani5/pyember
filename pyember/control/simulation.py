import sys
import os
import random
import copy
import yaml
from ..sim import MMCSim

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
            self.sim = MMCSim()
        elif self.sim_params['sim_type'] == 'KMC':
            self.sim = KMCSim()

        # Check compatibility of Config with Hamilton and Move objects
        if config is not None:
            _check_config_compatibility(config, hamilton, moves):


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


    def _check_config_compatibility(self, config, hamilton, moves):

        # Check Hamiltonian
        if config['latt_type'] != hamilton['latt_type']:
            raise ValueError("Config and Hamilton lattice types ara incompatible.")

        # Check Moves
        for key, move in moves.items():
            if config['latt_type'] != move['latt_type']:
                raise ValueError(f"Config and {key} move lattice types ara incompatible.")


    def setup(self, random_seed=42):
        """
        Initialize a selected model with parameters contained in a dict.
        """

        # read input configuration
        lat_type, box, xyz = read_cfg(self.incfg_file)

        # Initialize KMC system with appropriate lattice type
        self.kmc = Model(model_type, params, config)

        self.kmc.make_lattice(xyz, box)

        # read kMC parameters (reaction rates)
        rates = read_pars(self.param_file)

        # make event list (e.g., identify deposition sites)
        self.kmc.init_events(rates)

        # initialize random number generator
        random.seed(random_seed)


    def run(self, config):
        """
        Run simulation: call model to update configuration.
        In Metropolis MC, 'time' is measured in MC steps
        Simulation is stopped when final time is reached.

        Separate runs can be performed for equilibration and production.

        Parameters
        ----------
        config: dict or Config object
            initial configuration
        """


        # initial values
        t = t_print = t_save = t_measure = 0.0

        # print initial numbers
        print(t, self.kmc.nat)
        if (self.print_period < self.t_max):
            print('time, iteration, number of atoms')

        while t < self.t_max:

            t += 1.0
            it += 1

            # try move
            event = self.move(self.config)

            # energy difference
            du = self.hamilton(config, event)

            # accept move
            if du < 0:
                self.accept(self.config, event, self.hamilton)
            elif np.exp(-self.beta*du) > np.random.random():
                self.accept(self.config, event, self.hamilton)


            # perform runtime outputs
            if (t - t_print) > self.print_period:
                print(t, self.hamilton.get_energy_total(), 0.0)
                t_print = t

            if (t - t_save) > self.save_traj_period:
                io.write_xyz(self.config, 'mmc.xyz')
                t_save = t

            if (t - t_measure) > self.measure_period:
                t_measure = t


        print('End of simulation')


if __name__ == "__main__":

    # start simulation with parameters defined in an control file
    sim = Simulation(sys.argv[0])

    # initialize model and system configuration
    sim.setup()

    sim.run()

    sim.finalize()

