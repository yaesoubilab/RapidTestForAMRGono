from definitions import get_scenario_names, N_BREAKS_SENSITIVITY, N_BREAKS_SPECIFICITY
from model.Plots import plot_scenarios


plot_scenarios(scenario_names=get_scenario_names(n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
                                                 n_breaks_specificity=N_BREAKS_SPECIFICITY),
               fig_file_name='figures/Changing specificity.png')


