from SimulateScenarios import RAPID_TEST_COVERAGE
from definitions import get_scenario_names, N_BREAKS_SENSITIVITY, N_BREAKS_SPECIFICITY, MIN_SEN, MIN_SPE
from model.Plots import plot_scenarios

plot_scenarios(
    scenario_names=get_scenario_names(
        min_sensitivity=MIN_SEN,
        min_specificity=MIN_SPE,
        n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
        n_breaks_specificity=N_BREAKS_SPECIFICITY,
        rapid_test_coverage=RAPID_TEST_COVERAGE),
    fig_file_name='figures/Changing specificity.png')


