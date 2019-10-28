import sys
import os
import random

class KMCSim:
    """
    Class for simulation flow control of Metropolis Monte Carlo.
    """

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


    def _check_config(config_params):
        """Verifies if the control file configuration parameters are
        compatible with the config file.
        """

        config = read_xyz(sim_params['config']['file'])

        assert config['latt_type'] == config_params['type'], "Latt type in file does not match"
        assert config['pbc'] == config_params['pbc'], "PBC in file does not match"
        assert config['latt_box'] == config_params['latt_box'], "Latt box in file does not match"

        return config


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

