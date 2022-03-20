from definitions import get_scenario_names, N_BREAKS_SENSITIVITY, N_BREAKS_SPECIFICITY, MIN_SEN_SPE
from model.Plots import plot_scenarios


plot_scenarios(
    scenario_names=get_scenario_names(
        min_sensitivity=MIN_SEN_SPE,
        n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
        n_breaks_specificity=N_BREAKS_SPECIFICITY,
        n_breaks_rapid_test_coverage=1),
    fig_file_name='figures/Changing specificity.png')


