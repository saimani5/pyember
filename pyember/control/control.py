import sys
import os
import random

class Simulation:
    """
    Class for simulation flow control.
    """

    # recognized simulation flow control parameters, with optional defaults
    known_params = {
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

    def __init__(self, setup_info=None):
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
            self.sim_params = self.__read_control_dict(setup_info)
        elif isinstance(setup_info, str):
            self.sim_params = self.__read_control_file(setup_info)
        else:
            raise TypeError("Simulation setup info must be either a file name or dict")


    def __read_control_dict(self, setup_dict):
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
        for key, val in known_params.items():
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


    def __read_control_file(self, setup_file, directory='.'):
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


    def run(self):
        """
        Run simulation: call model to update configuration.
        In Metropolis MC, 'time' is measured in MC steps
        Simulation is stopped when final time is reached.

        Separate runs can be performed for equilibration and production.
        """

        # initial values
        t = t_print = t_save = t_measure = 0.0
        it = 0

        # print initial numbers
        if (self.print_period < self.t_max):
            print('time, iteration, number of atoms')
            print(t, it, self.kmc.nat)

        while t < self.t_max:

            t += self.kmc.advance_time()

            self.kmc.step()

            it += 1

            # perform runtime outputs
            if (t - t_print) > self.print_period:
                print(t, it, self.kmc.nat)
                t_print = t

            if (t - t_save) > self.save_traj_period:
                t_save = t

            if (t - t_measure) > self.measure_period:
                t_measure = t

        if (self.print_period < self.t_max):
            print('End of simulation')
            print(t, it, self.kmc.nat)


    def output(self):
        """
        Save the final state and statistics
        """

        xyz, box, grain = self.kmc.get_conf()
        write_cfg(self.outcfg_file, xyz, box, grain)


if __name__ == "__main__":

    # start simulation with parameters defined in an control file
    sim = Simulation(sys.argv[0])

    # initialize model and system configuration
    sim.setup()

    sim.run()

    sim.finalize()

