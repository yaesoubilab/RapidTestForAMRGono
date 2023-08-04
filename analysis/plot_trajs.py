import matplotlib

from definitions import get_scenario_name, get_traj_fig_name, SIM_DURATION
from model.plots import plot_trajectories

matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False


"""
To plot the simulated trajectories available under 'outputs/sim-with M-25yrs/trajectories'
"""


def plot_simulated_trajs(
        if_m_available=True, dict_test_characts=None, if_wider_priors=False, transmission_factor=1.0):

    # scenario name
    scenario_name = get_scenario_name(if_m_available=if_m_available,
                                      sim_duration=SIM_DURATION,
                                      calibration_seed=None,
                                      if_wider_priors=if_wider_priors)

    dir_of_traj_files = 'outputs/sim-{}/trajectories'.format(scenario_name)
    dir_of_traj_figs = 'figures/sim-{}/trajs'.format(scenario_name)

    figure_filename = get_traj_fig_name(
        if_m_available=if_m_available, dict_test_characts=dict_test_characts, transmission_factor=transmission_factor)

    # plot trajectories
    plot_trajectories(prev_multiplier=1,  # to show weeks on the x-axis of prevalence data
                      incd_multiplier=1,  # to show weeks on the x-axis of incidence data
                      obs_prev_multiplier=1,
                      obs_incd_multiplier=1,
                      dir_of_traj_files=dir_of_traj_files,
                      dir_of_traj_figs=dir_of_traj_figs,
                      filename=figure_filename)


# base (M is available)
plot_simulated_trajs(if_m_available=True)
