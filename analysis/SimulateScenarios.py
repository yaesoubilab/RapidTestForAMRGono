import apace.Calibration as calib
import model.Plots as P
from apace.ScenarioSimulation import ScenarioSimulator
from definitions import get_scenario_names, get_list_sens_spec_coverage, \
    N_BREAKS_SENSITIVITY, N_BREAKS_SPECIFICITY, MIN_SEN, MIN_SPE
from model.Model import build_model
from model.ModelSettings import GonoSettings

RAPID_TEST_COVERAGE = 0.75
IF_M_AVAILABLE_FOR_FIRST_TX = True
N_OF_SIMS = 160
RUN_IN_PARALLEL = True


def simulate_scenarios():

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=IF_M_AVAILABLE_FOR_FIRST_TX)
    sets.exportTrajectories = False

    # names of the scenarios to evaluate
    scenario_names = get_scenario_names(
        min_sensitivity=MIN_SEN,
        min_specificity=MIN_SPE,
        n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
        n_breaks_specificity=N_BREAKS_SPECIFICITY,
        rapid_test_coverage=RAPID_TEST_COVERAGE)

    # variable names (these correspond to the arguments of update_settings function of ModelSettings)
    var_names = ['sensitivity', 'specificity', 'rapid test coverage']

    # variable values
    # rows correspond to scenario names defined above, and columns correspond to variable names defined above
    # [0.0, 1.0, 0.0]  # status quo (no rapid test)
    scenario_definitions = [[0.0, 1.0, 0.0]] + get_list_sens_spec_coverage(
        min_sen=MIN_SEN,
        min_spe=MIN_SPE,
        n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
        n_breaks_specificity=N_BREAKS_SPECIFICITY,
        rapid_test_coverage=RAPID_TEST_COVERAGE)

    # get the seeds and probability weights
    seeds, lns, weights = calib.get_seeds_lnl_probs('outputs/calibration/calibration_summary.csv')

    scenario_sim = ScenarioSimulator(model_settings=sets,
                                     scenario_names=scenario_names,
                                     variable_names=var_names,
                                     scenario_definitions=scenario_definitions)

    scenario_sim.simulate(function_to_populate_model=build_model,
                          num_of_sims=N_OF_SIMS,
                          seeds=seeds, weights=weights, sample_seeds_by_weights=False,
                          if_run_in_parallel=RUN_IN_PARALLEL,
                          print_summary_stats=True)

    # export results of the scenario analysis
    scenario_sim.export_results()

    if IF_M_AVAILABLE_FOR_FIRST_TX:
        fig_file_name = 'figures/Changing specificity-with-M.png'
    else:
        fig_file_name = 'figures/Changing specificity-no-M.png'

    # plot the CEA figure and other analyses
    P.plot_scenarios(scenario_names=scenario_names,
                     csv_file_name=sets.folderToSaveScenarioAnalysis + '/simulated_scenarios.csv',
                     fig_file_name=fig_file_name)


if __name__ == "__main__":
    simulate_scenarios()
