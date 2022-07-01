import apacepy.calibration as calib

from model.model_settings import GonoSettings
from model.model_structure import build_model
from model.support import estimate_parameters, simulate_calibrated_model

"""
To calibrate the model under two scenarios: 
    1) when M is available for first-line therapy 
    2) when M is not available for first-line therapy

The results will be saved under analysis/outputs/with-M or analysis/outputs/no-M
"""

RUN_IN_PARALLEL = True
N_OF_CALIBRATION_ITERATIONS = 16*100    # total number of trajectories to simulate as part of calibration
N_OF_TRAJS_TO_USE_FOR_SIMULATION = 16*10   # number of trajectories with the highest likelihood to keep
N_OF_RESAMPLES_FOR_PARAM_ESTIMATION = 16*10  # number of parameter values to resample for parameter estimation


def calibrate(if_m_available):
    """ calibrate and simulate the calibrated model
    :param if_m_available: (bool) if M is available for first-line therapy
    """

    # get model settings
    sets = GonoSettings(if_calibrating=True, collect_traj_of_comparts=False,
                        if_m_available_for_1st_tx=if_m_available)
    # calibrate under the status quo scenario (no rapid test)
    sets.update_settings(sens=0, spec=1, prob_rapid_test=0)

    # --------- calibration ----------
    # calibrate the model
    calibration = calib.CalibrationWithRandomSampling(model_settings=sets, max_tries=200)

    calibration.run(
        function_to_populate_model=build_model,
        num_of_iterations=N_OF_CALIBRATION_ITERATIONS,
        if_run_in_parallel=RUN_IN_PARALLEL)

    # run time
    print('Run time: {:.2f} seconds for {} trajectories.'.format(calibration.runTime, N_OF_CALIBRATION_ITERATIONS))

    # save calibration results
    calibration.save_results()

    # estimate parameters and draw histograms
    estimate_parameters(n_of_resamples=N_OF_RESAMPLES_FOR_PARAM_ESTIMATION,
                        calibration_summary_file=sets.folderToSaveCalibrationResults+'/calibration_summary.csv',
                        if_m_available_for_1st_tx=if_m_available,
                        calibration_folder=sets.folderToSaveCalibrationResults)

    # simulate the calibrated model
    sets = GonoSettings(if_calibrating=False, collect_traj_of_comparts=True,
                        if_m_available_for_1st_tx=if_m_available)
    sets.update_settings(sens=0, spec=1, prob_rapid_test=0)

    if if_m_available:
        figure_name = 'Calibrated with M.png'
    else:
        figure_name = 'Calibrated no M.png'

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sample_seeds_by_weights=False,
                              settings=sets,
                              figure_filename=figure_name)


if __name__ == "__main__":

    calibrate(if_m_available=True)
    calibrate(if_m_available=False)
