import warnings

import apacepy.calibration as calib
from analyze_and_plot_scenarios import export_summary_and_plots_for_varying_coverage, plot_by_sens_spec
from apacepy.scenario_simulation import ScenarioSimulator
from definitions import get_sens_analysis_names_and_definitions, SIM_DURATION, TRANSMISSION_FACTOR_VALUES
from model.model_settings import GonoSettings
from model.model_structure import build_model

warnings.filterwarnings("ignore")

"""
To simulate and plot the impact of rapid tests with different characteristics
The results will be saved under outputs/(with or no)-M/scenarios
"""

N_OF_SIMS = 4 #250
RUN_IN_PARALLEL = True


def simulate_scenarios(if_m_available_for_1st_tx, simulation_duration,
                       vary_sens_spec=False, vary_transm_factor=False,
                       calibration_seed=None, if_wider_prior=False):
    """
    :param if_m_available_for_1st_tx:
    :param simulation_duration:
    :param vary_sens_spec: (bool) set true if the analysis should include varying
        the sensitivity and specificity of the test
    :param vary_transm_factor: (bool) set true if the analysis should include varying
        the transmission factor
    :param calibration_seed:
    :param if_wider_prior:
    :return:
    """

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available_for_1st_tx,
                        sim_duration=simulation_duration,
                        calibration_seed=calibration_seed,
                        if_wider_prior=if_wider_prior)
    sets.exportTrajectories = False

    # variable names (these correspond to the arguments of update_settings function of ModelSettings)
    var_names = ['CIP-sens', 'CIP-spec', 'TET-sens', 'TET-spec', 'rapid test coverage', 'transmission factor']

    # names of the scenarios to evaluate
    scenario_names, scenario_definitions = get_sens_analysis_names_and_definitions(
        vary_sens_spec=vary_sens_spec, vary_transm_factor=vary_transm_factor)

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
    if not vary_transm_factor:
        export_summary_and_plots_for_varying_coverage(
            if_m_available=if_m_available_for_1st_tx,
            simulation_duration=simulation_duration,
            calibration_seed=calibration_seed,
            trans_factor=1.0)
    else:
        for f in TRANSMISSION_FACTOR_VALUES:
            export_summary_and_plots_for_varying_coverage(
                if_m_available=if_m_available_for_1st_tx,
                simulation_duration=simulation_duration,
                calibration_seed=calibration_seed,
                trans_factor=f)

    if vary_sens_spec:
        plot_by_sens_spec(
            if_m_available=if_m_available_for_1st_tx,
            simulation_duration=simulation_duration,
            calibration_seed=calibration_seed)


if __name__ == "__main__":

    # print('\n*** M is available for 1st Tx ***')
    # simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=SIM_DURATION,
    #                    include_sens_analysis_on_sens_spec=True)

    print('\n*** M is available for 1st Tx with simulation duration of 35 years ***')
    simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=35)

    print('\n*** M is available for 1st Tx with varying transmission factor ***')
    simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=SIM_DURATION,
                       vary_transm_factor=True)

    # print('\n*** M is available for 1st Tx with a new initial calibration seed ***')
    # simulate_scenarios(if_m_available_for_1st_tx=True, simulation_duration=SIM_DURATION,
    #                    calibration_seed=1)
    #
    # print('\n*** M is not available for 1st Tx***')
    # simulate_scenarios(if_m_available_for_1st_tx=False, simulation_duration=SIM_DURATION)
    #
    # print('\n*** M is not available for 1st Tx with simulation duration of 35 years ***')
    # simulate_scenarios(if_m_available_for_1st_tx=False, simulation_duration=35)
