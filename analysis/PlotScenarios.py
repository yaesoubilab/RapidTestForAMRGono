from SimulateScenarios import RAPID_TEST_COVERAGE
from definitions import get_scenario_names, N_BREAKS_SENSITIVITY, N_BREAKS_SPECIFICITY, MIN_SEN, MIN_SPE
from model.Plots import plot_scenarios


IF_M_AVAILABLE_FOR_FIRST_TX = True

if IF_M_AVAILABLE_FOR_FIRST_TX:
    fig_file_name = 'figures/Changing specificity-with-M.png'
    csv_file_name = 'outputs-with-M/scenarios/simulated_scenarios.csv'
else:
    fig_file_name = 'figures/Changing specificity-no-M.png'
    csv_file_name = 'outputs-no-M/scenarios/simulated_scenarios.csv'


plot_scenarios(
    scenario_names=get_scenario_names(
        min_sensitivity=MIN_SEN,
        min_specificity=MIN_SPE,
        n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
        n_breaks_specificity=N_BREAKS_SPECIFICITY,
        rapid_test_coverage=RAPID_TEST_COVERAGE),
    csv_file_name=csv_file_name,
    fig_file_name=fig_file_name)


