import SimPy.ParameterEstimation as P
import apace.Calibration as calib
from apace.MultiEpidemics import MultiEpidemics
from definitions import ROOT_DIR
from model.Model import build_model
from model.ModelSettings import GonoSettings
from model.Plots import plot_trajectories


def simulate_multi_trajectories(n, seeds=None, weights=None, sample_seeds_by_weights=True,
                                if_run_in_parallel=True, figure_filename='traj.png',
                                settings=None):
    """
    :param n: (int) number of trajectories to simulate
    :param seeds: (list) of seeds
    :param weights: (list) probability weights over seeds
    :param sample_seeds_by_weights: (bool) set to False to only use seeds with positive weights
    :param if_run_in_parallel: (bool if run in parallel
    :param figure_filename: (string) filename to save the figures as
    :param settings: (GonoSettings) model settings
    :return:
    """

    # get model settings
    sets = GonoSettings() if settings is None else settings

    # build multiple epidemics
    multi_model = MultiEpidemics(model_settings=sets)

    multi_model.simulate(function_to_populate_model=build_model,
                         n=n,
                         seeds=seeds,
                         weights=weights,
                         sample_seeds_by_weights=sample_seeds_by_weights,
                         if_run_in_parallel=if_run_in_parallel)

    # save ids, seeds, runtime,
    multi_model.save_summary()

    # get summary statistics of runtime,
    multi_model.print_summary_stats()

    if sets.ifMAvailableFor1stTx:
        dir_of_trajs = 'outputs-with-M/trajectories'
    else:
        dir_of_trajs = 'outputs-no-M/trajectories'

    # plot trajectories
    plot_trajectories(prev_multiplier=1,  # to show weeks on the x-axis of prevalence data
                      incd_multiplier=1,  # to show weeks on the x-axis of incidence data
                      obs_prev_multiplier=1,
                      obs_incd_multiplier=1,
                      dir_of_trajs=dir_of_trajs,
                      filename=figure_filename)


def simulate_calibrated_model(n_of_sims, sample_seeds_by_weights=True,
                              if_run_in_parallel=True, settings=None,
                              figure_filename='Calibrated.png'):
    """
    simulates the calibrated model
    :param n_of_sims: (int) number of trajectories to simulate
    :param sample_seeds_by_weights: (bool)
    :param if_run_in_parallel: (bool if run in parallel
    :param settings: (GonoSettings) model settings
    :param figure_filename: (string) figure name
    """

    # ------- simulate the calibrated model ----------
    # get the seeds and probability weights
    seeds, ln, weights = calib.get_seeds_lnl_probs('outputs/calibration/calibration_summary.csv')

    # simulate the calibrated model
    simulate_multi_trajectories(n=n_of_sims,
                                seeds=seeds,
                                weights=weights,
                                sample_seeds_by_weights=sample_seeds_by_weights,
                                if_run_in_parallel=if_run_in_parallel,
                                figure_filename=figure_filename,
                                settings=settings)


def estimate_parameters(n_of_resamples):
    """
    :param n_of_resamples:(int) number of parameter values to resamples
    """
    # ---------- parameter estimation -----------
    # calculate posterior distributions and plot figures
    estimator = P.ParameterAnalyzer()
    estimator.resample_param_values(
        csvfile_param_values_and_weights='outputs/calibration/calibration_summary.csv',
        n=n_of_resamples,
        weight_col=3,
        sample_by_weight=False,
        csvfile_resampled_params='outputs/calibration/resampled_parameter_values.csv',
        seed=0)

    param_list_for_table = [
        'Transmission parameter',
        'Time until natural recovery',
        'Time until screened',
        'Time until seeking treatment (symptomatic)',
        'Time until seeking retreatment (symptomatic)',
        'Prob symptomatic',

        'Exponent for the prob of resistance by antibiotics-0',
        'Exponent for the prob of resistance by antibiotics-1',
        'Exponent for the prob of resistance by antibiotics-2'] \
        + ['Relative infectivity by infectivity profile-{}'.format(i) for i in range(8)] \
        + ['Initial prevalence',
           'Initial % I by symptom states-0'] \
        + ['Initial % I by resistance profile-{}'.format(i) for i in range(8)]

    param_list_for_figure = [
        'Transmission parameter',
        'Time until natural recovery',
        'Time until screened',
        'Relative infectivity by infectivity profile-1',
        'Relative infectivity by infectivity profile-2',
        'Relative infectivity by infectivity profile-3',
        'Prob symptomatic',
        'Exponent for the prob of resistance by antibiotics-0',
        'Exponent for the prob of resistance by antibiotics-1',
        'Exponent for the prob of resistance by antibiotics-2'
    ]
    # print('\nPosterior distributions:')
    # estimator.print_means_and_intervals(param_names=param_list_for_table)
    estimator.export_means_and_intervals(poster_file='outputs/calibration/posteriors.csv',
                                         param_names=param_list_for_table,
                                         prior_info_csv_file=ROOT_DIR+'/model/data/priors.csv')
    estimator.plot_pairwise(fig_filename='figures/posterior_figure.png',
                            par_names=param_list_for_figure,
                            prior_info_csv_file=ROOT_DIR + '/model/data/priors.csv')
