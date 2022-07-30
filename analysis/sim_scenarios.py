import warnings

import apacepy.calibration as calib
from analyze_and_plot_scenarios import export_performance_summary_and_plots
from apacepy.scenario_simulation import ScenarioSimulator
from definitions import get_sens_analysis_names_and_definitions, SIM_DURATION
from model.model_settings import GonoSettings
from model.model_structure import build_model

warnings.filterwarnings("ignore")

"""
To simulate and plot the impact of rapid tests with different characteristics
The results will be saved under outputs/(with or no)-M/scenarios
"""

N_OF_SIMS = 16
RUN_IN_PARALLEL = True


def simulate_scenarios(if_m_available_for_1st_tx, simulation_duration,
                       include_sens_analysis_on_sens_spec=False, calibration_seed=None):

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available_for_1st_tx,
                        sim_duration=simulation_duration,
                        calibration_seed=calibration_seed)
    sets.exportTrajectories = False

    # variable names (these correspond to the arguments of update_settings function of ModelSettings)
    var_names = ['CIP-sens', 'CIP-spec', 'TET-sens', 'TET-spec', 'rapid test coverage']

    # names of the scenarios to evaluate
    scenario_names, scenario_definitions = get_sens_analysis_names_and_definitions(
        include_sens_analysis_on_sens_spec=include_sens_analysis_on_sens_spec)

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
    simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=SIM_DURATION)

    print('\n*** M is available for 1st Tx with simulation duration of 30 years ***')
    simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=35)

    print('\n*** M is available for 1st Tx with a new initial calibration seed ***')
    simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=SIM_DURATION,
                       calibration_seed=1)

    print('\n*** M is not available for 1st Tx***')
    simulate_scenarios(if_m_available_for_1st_tx=False, simulation_duration=SIM_DURATION)

    print('\n*** M is not available for 1st Tx with simulation duration of 30 years ***')
    simulate_scenarios(if_m_available_for_1st_tx=False, simulation_duration=35)
