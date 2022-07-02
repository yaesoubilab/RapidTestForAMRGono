from calibrate import N_OF_RESAMPLES_FOR_PARAM_ESTIMATION
from definitions import get_scenario_name
from model.support import estimate_parameters

"""
To estimate posterior distribution of model parameters after calibrating the model
and to plot the the posterior distribution and correlation between parameters 
"""


def estiamte_parameter(if_m_available, calibration_seed):

    # scenario name
    scenario_name = get_scenario_name(if_m_available=if_m_available, calibration_seed=calibration_seed)

    calibration_folder = 'outputs/{}/calibration'.format(scenario_name_without_sim_duration)

    estimate_parameters(n_of_resamples=N_OF_RESAMPLES_FOR_PARAM_ESTIMATION,
                        calibration_summary_file=sets.folderToSaveCalibrationResults + '/calibration_summary.csv',
                        calibration_folder=sets.folderToSaveCalibrationResults,
                        figure_folder='figures/' + scenario_name
                        )

