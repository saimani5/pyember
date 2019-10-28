import numpy as np
from ..io import read_xyz, write_xyz
from ..interact import Heisenberg
from ..move import MMCMove

class MMCSim:
    """
    Class for simulation flow control of Metropolis Monte Carlo.
    """

    # recognized simulation flow control parameters, with optional defaults
    def __init__(self, sim_params):
        """
        Initializes simulation object either from an input file or a
        dictionary with appropriate parameters.

        Parameters
        ----------
        sim_params: dict
            contains information for setting up a simulation flow control
        """

        self.mmc_params = {}
        self.mmc_params['time_control'] = sim_params['time_control']
        self.mmc_params['moves'] = sim_params['moves']
        config = self._check_config(sim_params['config'])
        self.mmc_params['config'] = config
        self.config = config
        self.t_max = self.mmc_params['time_control']['total']
        self.print_period = self.mmc_params['time_control']['print']
        self.save_traj_period = self.mmc_params['time_control']['save']
        self.measure_period = self.mmc_params['time_control']['measure']

        # Supported Hamiltonians
        hamilton = {'heisenberg': Heisenberg}

        # Set up Hamiltonian and check its compatibility with config
        ham_type = sim_params['hamilton']['type']
        ham_params = sim_params['hamilton']['params']
        ham = hamilton[ham_type](config, ham_params)
        self.mmc_params['hamilton'] = ham
        self.hamilton = ham
        self.du = ham.get_energy_diff_i
        self.get_energy_total = ham.get_energy_total

        # Set up moves
        self.mmc_params['moves'] = MMCMove(sim_params['moves'], config)
        self.move = self.mmc_params['moves'].move
        self.accept = self.mmc_params['moves'].accept

        # initialize random number generator
        np.random.seed(sim_params['random_seed'])



    def _check_config(self, config_params):
        """Verifies if the control file configuration parameters are
        compatible with the config file.
        """
        config = read_xyz(config_params['file'])

        assert config['latt_type'] == config_params['type'], "Latt type in file does not match"
        assert config['pbc'] == config_params['pbc'], f"PBC in file does not match {config['pbc']} vs {config_params['pbc']}"
        assert np.array_equal(config['latt_box'], config_params['latt_box']), "Latt box in file does not match"

        return config


    def run(self, config=None):
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

        if config is not None:
            self.config = config

        # initial values
        t = t_print = t_save = t_measure = 0.0

        # initial energy and magnetization statistics
        tot_ene = self.hamilton.get_energy_total(self.config)
        tot_mag, tsx, tsy, tsz = self.hamilton.get_magnetization(self.config)
        print('time, total_energy, |M|, Mx, My, Mz')
        print(t, tot_ene, round(tot_mag), round(tsx), round(tsy), round(tsz))

        while t < self.t_max:

            t += 1.0

            # try move
            event = self.move(self.config)

            # energy difference
            beta_du = self.du(self.config, event)
            #print('du', beta_du)

            # accept move
            if beta_du < 0:
                self.accept(self.config, event, self.hamilton)
            elif np.exp(-beta_du) > np.random.random():
                self.accept(self.config, event, self.hamilton)


            # perform runtime outputs
            if (t - t_print) > self.print_period:
                tot_ene = self.hamilton.get_energy_total(self.config)
                tot_mag, tsx, tsy, tsz = self.hamilton.get_magnetization(self.config)
                print(t, tot_ene, round(tot_mag), round(tsx), round(tsy), round(tsz))
                t_print = t

            if (t - t_save) > self.save_traj_period:
                write_xyz(self.config, 'mmc.xyz')
                t_save = t

            if (t - t_measure) > self.measure_period:
                t_measure = t

        print('End of simulation')
