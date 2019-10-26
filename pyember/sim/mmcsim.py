import sys
import os
import random

class MMCSim(Sim):
    """
    Class for simulation flow control of Metropolis Monte Carlo.
    """

    # recognized simulation flow control parameters, with optional defaults
    self.known_params = {
            't_max':100.0,
            'print_period':1,
            'save_period':100,
            'measure_period':100,
            'model_params_file':'model.params',
            'incfg_file':'input.xyz',
            'outcfg_file':None,
            'traj_file':None,
            'stats_file':None
        }

    def __init__(self, setup_info=None, hamilton=None, moves=None, config=None):
        """
        Initializes simulation object either from an input file or a
        dictionary with appropriate parameters.

        Parameters
        ----------
        setup_info: str or dict
            contains information for setting up a simulation flow control
            variables, which are then stored in a sim_parms dict.
        """

        if isinstance(setup_info, dict):
            self.sim_params = self._read_control_dict(setup_info)
        elif isinstance(setup_info, str):
            self.sim_params = self._read_control_file(setup_info)
        else:
            raise TypeError("Simulation setup info must be either a file name or dict")

        if hamilton is None:
            raise ValueError("Missing Hamiltonian object")

        if moves is None or moves == {}:
            raise ValueError("Missing Move objects")

        # Check compatibility of Config with Hamilton and Move objects
        if config is not None:
            _check_config_compatibility(config, hamilton, moves):


    def _read_control_dict(self, setup_dict):
        """
        Verifies that the supplied dictionary contains the necessary parameters
        and no unknown parameters.

        Parameters
        ----------
        setup_dict: dict
            Supplied dict of parameter values

        Returns
        -------
        param_dict: dict
            Validated dict of parameter values
        """

        # check if all necessary parameters are present
        for key, val in self.known_params.items():
            if val and (key not in setup_dict.keys()):
                raise ValueError(f"{key} is required, but not present in the supplied parameters")

        param_dict = {}

        # check if all supplied parameters are meaningful
        for key in setup_dict:
            try:
                param_dict[key] = setup_dict[key]
            except KeyError as e:
                print(f"{key}: Unknown simulation parameter")
                raise e

        return param_dict


    def _read_control_file(self, setup_file, directory='.'):
        """
        Read control input file into a dict.
        The file is structured as key:value pairs in no particular order. The
        key is used in the parameter dict.

        Parameters
        ----------
        setup_file: str
            file containing information for setting up the simulation
        directory: str
            if supplied, it provides the default directory for input and
            output files, defaults to the current directory

        Returns
        -------
        param_dict: dict
            dict of parameters and their values
        """
        
        setup_dict = {}

        with open(setup_file, 'r') as f:
            for line in iter(f.readline, ''):
                key, value = line.split(':')
                key = key.strip()
                value = value.strip()
                setup_dict[key] = value

        # validate the parameter dict
        param_dict = self.__read_control_dict(setup_dict)

        return param_dict

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

    self.known_params = {
            't_max':100.0,
            'print_period':1,
            'save_period':100,
            'measure_period':100,
        }

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

