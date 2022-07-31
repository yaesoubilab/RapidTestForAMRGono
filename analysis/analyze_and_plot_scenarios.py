from definitions import COVERAGE_VALUES, SIM_DURATION
from model.plots import plot_sa_for_a_specific_coverage, plot_sa_for_varying_coverage, plot_scenario_sa
from model.scenario_and_sensitivity_analyses import get_scenarios_csv_filename_and_fig_filename, \
    export_performance_of_scenarios

"""
To plot the cost-effectiveness plan visualizing the performance of rapid tests for 
different test characteristics. The results will be save under 
    analysis/figures/SA and 
    analysis/outputs/(with or no-M)/scenarios/performance_summary.csv
"""


X_RANGE_WITH_M = [-0.1, 6.1]
Y_RANGE_WITH_M = [-250, 1750]
X_RANGE_NO_M = [-0.1, 6.1]
Y_RANGE_NO_M = [-4000, 500]


def plot_sa_for_a_test_coverage(
        if_m_available_for_1st_tx, test_coverage, sim_duration=None, calibration_seed=None):
    """
    plots a cost-effectiveness plan visualizing the performance of a rapid test for
        the specified scenario for the availability of M and test coverage
    :param if_m_available_for_1st_tx: (bool) if M is available for the first-line therapy
    :param test_coverage: (float) specified test coverage
    :param sim_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    """

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available_for_1st_tx,
        test_coverage=test_coverage,
        sim_duration=sim_duration,
        calibration_seed=calibration_seed)

    if if_m_available_for_1st_tx:
        x_range = X_RANGE_WITH_M
        y_range = Y_RANGE_WITH_M
    else:
        x_range = X_RANGE_NO_M
        y_range = Y_RANGE_NO_M

    plot_sa_for_a_specific_coverage(
        csv_file_name=csv_file_name,
        fig_file_name=fig_file_name,
        test_coverage=test_coverage,
        x_range=x_range,
        y_range=y_range)


def plot_scenarios_sa(if_m_available_for_1st_tx, include_sens_analysis_on_sens_spec=False,
                      sim_duration=None, calibration_seed=None):
    """
    plots a 3-panel figure of cost-effectiveness planes visualizing the performance of the rapid test for
        the specified scenario for the availability of M
    :param if_m_available_for_1st_tx: (bool) if M is available for the first-line therapy
    :param include_sens_analysis_on_sens_spec: (bool) if include analyses to vary sensitivity and specificity
        of CIP and TET
    :param sim_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    """

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available_for_1st_tx,
        sim_duration=sim_duration,
        calibration_seed=calibration_seed
    )

    if if_m_available_for_1st_tx:
        x_range = X_RANGE_WITH_M
        y_range = Y_RANGE_WITH_M
    else:
        x_range = X_RANGE_NO_M
        y_range = Y_RANGE_NO_M

    plot_scenario_sa(csv_file_name=csv_file_name, fig_file_name=fig_file_name,
                     x_range=x_range, y_range=y_range,
                     l_b_r_t=(0.1, 0.15, 0.97, 0.9),
                     fig_size=(8, 3.5))


def export_summary_and_plots_for_varying_coverage(if_m_available,
                                                  simulation_duration, calibration_seed=None):
    """
    exports performance summary and plot the cost-effectiveness figure for varying coverage level under
        this scenario of the availability for M
    :param if_m_available: (bool) if M is available for the first-line therapy
    :param simulation_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    """

    # export performance of different scenarios of test characteristics
    export_performance_of_scenarios(
        if_m_available_for_1st_tx=if_m_available,
        simulation_duration=simulation_duration,
        calibration_seed=calibration_seed,
        coverage_values=COVERAGE_VALUES)

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available,
        sim_duration=simulation_duration,
        calibration_seed=calibration_seed)

    if if_m_available:
        x_range = X_RANGE_WITH_M
        y_range = Y_RANGE_WITH_M
    else:
        x_range = X_RANGE_NO_M
        y_range = Y_RANGE_NO_M

    plot_sa_for_varying_coverage(
        csv_file_name=csv_file_name,
        sim_duration=simulation_duration,
        fig_file_name=fig_file_name,
        x_range=x_range,
        y_range=y_range)


def plot_by_sens_spec_coverage(if_m_available, simulation_duration=None, calibration_seed=None):
    """
    plot cost-effectiveness plans for this scenario for the availability of M
    :param if_m_available: (bool) if M is available for the first-line therapy
    :param simulation_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    """

    #  plots a 3-panel figure of cost-effectiveness planes visualizing the performance of the rapid test for
    #  the specified scenario for the availability of M
    # plot_scenarios_sa(
    #     if_m_available_for_1st_tx=if_m_available,
    #     include_sens_analysis_on_sens_spec=include_sens_analysis_on_sens_spec,
    #     sim_duration=simulation_duration,
    #     calibration_seed=calibration_seed
    # )

    for c in COVERAGE_VALUES:
        # plot a cost-effectiveness plan visualizing the performance of the rapid test for
        # the specified scenario for the availability of M and test coverage
        plot_sa_for_a_test_coverage(
            if_m_available_for_1st_tx=if_m_available,
            sim_duration=simulation_duration,
            calibration_seed=calibration_seed,
            test_coverage=c)


if __name__ == "__main__":

    export_summary_and_plots_for_varying_coverage(if_m_available=True, simulation_duration=SIM_DURATION)


    # # for scenarios where drug M is and is not available for 1st-line therapy
    # for m_available in (True, False):
    #     export_performance_summary_and_plots(if_m_available=m_available)

