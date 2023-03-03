from calibrate import N_OF_RESAMPLES_FOR_PARAM_ESTIMATION
from definitions import get_scenario_name, ROOT_DIR
from model.support import estimate_parameters

"""
To estimate posterior distribution of model parameters after calibrating the model
and to plot the the posterior distribution and correlation between parameters 
"""


def estimate_parameters_of_a_scenario(if_m_available, calibration_seed):

    # scenario name
    scenario_name = get_scenario_name(
        if_m_available=if_m_available, calibration_seed=calibration_seed)

    calibration_folder = ROOT_DIR+'/analysis/outputs/{}/calibration'.format(scenario_name)

    estimate_parameters(n_of_resamples=N_OF_RESAMPLES_FOR_PARAM_ESTIMATION,
                        calibration_summary_file=calibration_folder + '/calibration_summary.csv',
                        calibration_folder=calibration_folder,
                        figure_folder='figures/' + scenario_name
                        )


estimate_parameters_of_a_scenario(if_m_available=True, calibration_seed=None)