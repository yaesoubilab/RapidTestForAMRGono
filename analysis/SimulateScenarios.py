import warnings

import apace.Calibration as calib
import model.Plots as P
from analysis.PlotScenarios import X_RANGE_WITH_M, Y_RANGE_WITH_M
from apace.ScenarioSimulation import ScenarioSimulator
from definitions import get_scenario_names, get_list_sens_spec_coverage, COVERAGE_VALUES
from model.Model import build_model
from model.ModelSettings import GonoSettings

warnings.filterwarnings("ignore")


IF_M_AVAILABLE_FOR_FIRST_TX = True
N_OF_SIMS = 160
RUN_IN_PARALLEL = True


def simulate_scenarios(if_m_available_for_1st_tx):

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available_for_1st_tx)
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
    seeds, lns, weights = calib.get_seeds_lnl_probs('outputs/calibration/calibration_summary.csv')

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

    # plot the CEA figure and other analyses
    for c in COVERAGE_VALUES:

        if if_m_available_for_1st_tx:
            fig_file_name = 'figures/SA-with M-coverage {:.2f}.png'.format(c)
        else:
            fig_file_name = 'figures/SA-no M-coverage {:.2f}.png'.format(c)

        P.plot_scenarios(
            csv_file_name=sets.folderToSaveScenarioAnalysis + '/simulated_scenarios.csv',
            fig_file_name=fig_file_name,
            test_coverage=c,
            x_range=X_RANGE_WITH_M,
            y_range=Y_RANGE_WITH_M)


if __name__ == "__main__":

    print('\n*** M is available for 1st Tx ***')
    simulate_scenarios(if_m_available_for_1st_tx=True)

    print('\n*** M is not available for 1st Tx***')
    simulate_scenarios(if_m_available_for_1st_tx=False)
