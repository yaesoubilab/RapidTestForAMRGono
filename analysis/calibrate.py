import apacepy.calibration as calib
from deampy.in_out_functions import TextFile

from definitions import get_scenario_name, ROOT_DIR, SIM_DURATION
from model.model_settings import GonoSettings
from model.model_structure import build_model
from model.support import estimate_parameters, simulate_calibrated_model

"""
To calibrate the model under scenarios: 
    1) when M is available for first-line therapy 
    2) when a different seed is used
    3) when winder prior distributions are used 

"""

RUN_IN_PARALLEL = True
N_OF_CALIBRATION_ITERATIONS = 16*100    # total number of trajectories to simulate as part of calibration
N_OF_TRAJS_TO_USE_FOR_SIMULATION = 200   # number of trajectories with the highest likelihood to keep
N_OF_RESAMPLES_FOR_PARAM_ESTIMATION = 200  # number of parameter values to resample for parameter estimation


def calibrate(if_m_available, calibration_seed, if_wider_priors=False):
    """ calibrate and simulate the calibrated model
    :param if_m_available: (bool) if M is available for first-line therapy
    :param calibration_seed: (int or None) calibration seed (for sensitivity analysis)
    :param if_wider_priors: (bool) set to True if the winder prior distributions should be used
    """

    # scenario name
    scenario_name = get_scenario_name(
        if_m_available=if_m_available, calibration_seed=calibration_seed,
        sim_duration=None, if_wider_priors=if_wider_priors)

    # get model settings
    sets = GonoSettings(if_calibrating=True, collect_traj_of_comparts=False,
                        if_m_available_for_1st_tx=if_m_available,
                        calibration_seed=calibration_seed,
                        if_wider_priors=if_wider_priors)
    # calibrate under the status quo scenario (no rapid test)
    sets.update_settings(cip_sens=0, cip_spec=1, tet_sens=0, tet_spec=1, prob_rapid_test=0)

    # --------- calibration ----------
    # calibrate the model
    calibration = calib.CalibrationWithRandomSampling(model_settings=sets, max_tries=200)

    calibration.run(
        function_to_populate_model=build_model,
        num_of_iterations=N_OF_CALIBRATION_ITERATIONS,
        initial_seed=calibration_seed,
        if_run_in_parallel=RUN_IN_PARALLEL)

    # run time
    print('Run time: {:.2f} seconds for {} trajectories.'.format(calibration.runTime, N_OF_CALIBRATION_ITERATIONS))

    # store summary of calibration
    file = TextFile(filename=ROOT_DIR+'/analysis/outputs/calib-{}/calibration/stats.txt'.format(scenario_name))
    file.write('Number of calibration iterations: {}\n'.format(N_OF_CALIBRATION_ITERATIONS))
    file.write('Number of trajectories discarded: {}\n'.format(calibration.nTrajsDiscarded))
    file.write('Calibration duration (seconds): {}\n'.format(round(calibration.runTime, 1)))
    file.write('Number of trajectories with non-zero probability: {}\n'.format(calibration.nTrajsWithNonZeroProb))
    file.close()

    # save calibration results
    calibration.save_results()

    # estimate parameters and draw histograms
    estimate_parameters(n_of_resamples=N_OF_RESAMPLES_FOR_PARAM_ESTIMATION,
                        calibration_summary_file=sets.folderToSaveCalibrationResults+'/calibration_summary.csv',
                        calibration_folder=sets.folderToSaveCalibrationResults,
                        figure_folder='figures/calib-'+scenario_name)

    # simulate the calibrated model
    sets = GonoSettings(if_calibrating=False, collect_traj_of_comparts=True,
                        if_m_available_for_1st_tx=if_m_available,
                        calibration_seed=calibration_seed,
                        if_wider_priors=if_wider_priors)
    sets.update_settings(cip_sens=0, cip_spec=1, tet_sens=0, tet_spec=1, prob_rapid_test=0)

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sim_duration=SIM_DURATION,
                              calibration_seed=calibration_seed,
                              sample_seeds_by_weights=False,
                              settings=sets)


if __name__ == "__main__":

    # calibrate the model under the scenario when M is available for the first-line therapy
    # for two different seeds
    for calib_seed in [None, 1]:

        scen_name = get_scenario_name(
            if_m_available=True, calibration_seed=calib_seed, sim_duration=None)
        print("\nCalibrating scenario '{}':".format(scen_name))

        calibrate(if_m_available=True, calibration_seed=calib_seed)

    # calibrating the model using wider prior distributions
    scen_name = get_scenario_name(
        if_m_available=True, calibration_seed=None, sim_duration=None, if_wider_priors=True)
    print("\nCalibrating scenario '{}':".format(scen_name))

    calibrate(if_m_available=True, calibration_seed=None, if_wider_priors=True)
