import model.Plots as P
from apace.ScenarioSimulation import ScenarioSimulator
from definitions import get_scenario_names, get_list_sensitivity_specificity, \
    N_BREAKS_SENSITIVITY, N_BREAKS_SPECIFICITY
from model.Model import build_model
from model.ModelSettings import GonoSettings

N_OF_SIMS = 100
RUN_IN_PARALLEL = True


def simulate_scenarios():

    # get model settings
    sets = GonoSettings()
    sets.exportTrajectories = False

    # names of the scenarios to evaluate
    scenario_names = get_scenario_names(n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
                                        n_breaks_specificity=N_BREAKS_SPECIFICITY)

    # variable names (these correspond to the arguments of update_settings function of ModelSettings)
    var_names = ['sensitivity', 'specificity']

    # variable values
    # rows correspond to scenario names defined above, and columns correspond to variable names defined above
    scenario_definitions = get_list_sensitivity_specificity(n_breaks_sensitivity=N_BREAKS_SENSITIVITY,
                                                            n_breaks_specificity=N_BREAKS_SPECIFICITY)
    scenario_sim = ScenarioSimulator(model_settings=sets,
                                     scenario_names=scenario_names,
                                     variable_names=var_names,
                                     scenario_definitions=scenario_definitions)

    scenario_sim.simulate(function_to_populate_model=build_model,
                          num_of_sims=N_OF_SIMS,
                          if_run_in_parallel=RUN_IN_PARALLEL,
                          print_summary_stats=True)

    # export results of the scenario analysis
    scenario_sim.export_results()

    # plot the CEA figure and other analyses
    P.plot_scenarios(scenario_names=scenario_names,
                     fig_file_name='figures/Changing specificity.png')


if __name__ == "__main__":
    simulate_scenarios()
