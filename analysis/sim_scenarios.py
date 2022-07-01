import warnings

import apacepy.calibration as calib
from analyze_and_plot_scenarios import export_performance_summary_and_plots
from apacepy.scenario_simulation import ScenarioSimulator
from definitions import get_scenario_names, get_list_sens_spec_coverage
from model.model_settings import GonoSettings
from model.model_structure import build_model

warnings.filterwarnings("ignore")

"""
To simulate and plot the impact of rapid tests with different characteristics
The results will be saved under outputs/(with or no)-M/scenarios
"""

N_OF_SIMS = 2
RUN_IN_PARALLEL = False


def simulate_scenarios(if_m_available_for_1st_tx, simulation_duration=None, calibration_seed=None):

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available_for_1st_tx,
                        sim_duration=simulation_duration,
                        calibration_seed=calibration_seed)
    sets.exportTrajectories = False

    # names of the scenarios to evaluate
    scenario_names = get_scenario_names()

    # variable names (these correspond to the arguments of update_settings function of ModelSettings)
    var_names = ['sensitivity', 'specificity', 'rapid test coverage']

    # variable values
    # rows correspond to scenario names defined above, and columns correspond to variable names defined above
    # [0.0, 1.0, 0.0]  # status quo (no rapid test)
    scenario_definitions = [[0.0, 1.0, 0.0]] + get_list_sens_spec_coverage()

    # get the seeds and probability weights
    seeds, lns, weights = calib.get_seeds_lnl_probs(sets.folderToSaveCalibrationResults+'/calibration_summary.csv')

    scenario_sim = ScenarioSimulator(model_settings=sets,
                                     scenario_names=scenario_names,
                                     variable_names=var_names,
                                     scenario_definitions=scenario_definitions)

    scenario_sim.simulate(function_to_populate_model=build_model,
                          num_of_sims=N_OF_SIMS,
                          seeds=seeds, weights=weights, sample_seeds_by_weights=False,
                          if_run_in_parallel=RUN_IN_PARALLEL,
                          print_summary_stats=False)

    # export results of the scenario analysis
    scenario_sim.export_results()

    # export the summary of performance and cost-effectiveness plots
    export_performance_summary_and_plots(
        if_m_available=if_m_available_for_1st_tx,
        simulation_duration=simulation_duration,
        calibration_seed=calibration_seed
    )


if __name__ == "__main__":

    print('\n*** M is available for 1st Tx ***')
    simulate_scenarios(if_m_available_for_1st_tx=True)

    print('\n*** M is available for 1st Tx with simulation duration of 30 years ***')
    simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=30)

    print('\n*** M is available for 1st Tx with a new initial calibration seed ***')
    simulate_scenarios(if_m_available_for_1st_tx=True, calibration_seed=1)

    print('\n*** M is not available for 1st Tx***')
    simulate_scenarios(if_m_available_for_1st_tx=False)

    print('\n*** M is not available for 1st Tx with simulation duration of 30 years ***')
    simulate_scenarios(if_m_available_for_1st_tx=False, simulation_duration=30)
