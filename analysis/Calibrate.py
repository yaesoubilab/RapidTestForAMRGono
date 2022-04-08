import apace.Calibration as calib
from model import Model as M
from model.ModelSettings import GonoSettings
from model.Support import estimate_parameters, simulate_calibrated_model


RUN_IN_PARALLEL = True
N_OF_CALIBRATION_ITERATIONS = 160*10    # total number of trajectories to simulate as part of calibration
N_OF_TRAJS_TO_USE_FOR_SIMULATION = 16*10   # number of trajectories with the highest likelihood to keep
N_OF_RESAMPLES_FOR_PARAM_ESTIMATION = 16*10  # number of parameter values to resample for parameter estimation

if __name__ == "__main__":

    # get model settings
    sets = GonoSettings(if_calibrating=True, collect_traj_of_comparts=False)

    # --------- calibration ----------
    # calibrate the model
    calibration = calib.CalibrationWithRandomSampling(model_settings=sets, max_tries=200)

    calibration.run(
        function_to_populate_model=M.build_model,
        num_of_iterations=N_OF_CALIBRATION_ITERATIONS,
        if_run_in_parallel=RUN_IN_PARALLEL)

    # run time
    print('Run time: {:.2f} seconds for {} trajectories.'.format(calibration.runTime, N_OF_CALIBRATION_ITERATIONS))

    # save calibration results
    calibration.save_results(filename='outputs/calibration/calibration_summary.csv')

    # estimate parameters and draw histograms
    estimate_parameters(n_of_resamples=N_OF_RESAMPLES_FOR_PARAM_ESTIMATION)

    # simulate the calibrated model
    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sample_seeds_by_weights=False)
