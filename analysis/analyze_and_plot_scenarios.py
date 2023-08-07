import matplotlib

from definitions import COVERAGE_VALUES, TRIPLE_FIG_SIZE, SIM_DURATION
from model.plots import plot_sa_for_specific_ab, plot_sa_for_varying_coverage, plot_sa_for_specific_ab_and_coverage
from model.scenario_and_sensitivity_analyses import get_scenarios_csv_filename_and_fig_filename, \
    export_performance_of_scenarios

matplotlib.use('cairo')
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False


"""
To plot the cost-effectiveness plan visualizing the performance of rapid tests for 
different test characteristics. The results will be save under 
    analysis/figures/SA and 
    analysis/outputs/(with or no-M)/scenarios/performance_summary.csv
"""


X_RANGE_WITH_M = [-0.1, 6.1]
Y_RANGE_WITH_M = [-250, 2000]
X_RANGE_NO_M = [-0.1, 6.1]
Y_RANGE_NO_M = [-5000, 500]
X_RANGE_VARYING_F = [-0.1, 6.1]
Y_RANGE_VARYING_F = [-2000, 3000]


def plot_by_sens_spec_coverage(if_m_available_for_1st_tx, ab, include_sens_labels,
                               sim_duration=None, calibration_seed=None):
    """
    plots a 3-panel figure of cost-effectiveness planes visualizing the performance of the rapid test for
        the specified scenario for the availability of M
    :param if_m_available_for_1st_tx: (bool) if M is available for the first-line therapy
    :param ab: (string) 'CIP' or 'TET'
    :param include_sens_labels: (bool) set to True to show the sensitivity values on the figure
    :param sim_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    """

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available_for_1st_tx,
        ab=ab,
        sim_duration=sim_duration,
        calibration_seed=calibration_seed
    )

    if if_m_available_for_1st_tx:
        x_range = X_RANGE_WITH_M
        y_range = Y_RANGE_WITH_M
    else:
        x_range = X_RANGE_NO_M
        y_range = Y_RANGE_NO_M

    plot_sa_for_specific_ab(
        ab=ab, include_sens_labels=include_sens_labels,
        csv_file_name=csv_file_name, fig_file_name=fig_file_name,
        x_range=x_range, y_range=y_range,
        l_b_r_t=(0.1, 0.15, 0.97, 0.9), fig_size=TRIPLE_FIG_SIZE)


def export_summary_and_plots_for_varying_coverage(
        if_m_available, simulation_duration,
        calibration_seed=None, varying_trans_factor=False, if_wider_priors=False, interval='c'):
    """
    exports performance summary and plot the cost-effectiveness figure for varying coverage level under
        this scenario of the availability for M
    :param if_m_available: (bool) if M is available for the first-line therapy
    :param simulation_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    :param varying_trans_factor: (float) transmission factor
    :param if_wider_priors: (bool) if use the model with wider priors
    :param interval: (string) 'c' for confidence interval and 'p' for prediction interval
    """

    # export performance of different scenarios of test characteristics
    export_performance_of_scenarios(
        if_m_available_for_1st_tx=if_m_available,
        simulation_duration=simulation_duration,
        calibration_seed=calibration_seed,
        coverage_values=COVERAGE_VALUES,
        varying_trans_factor=varying_trans_factor,
        if_wider_priors=if_wider_priors)

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available,
        sim_duration=simulation_duration,
        calibration_seed=calibration_seed,
        varying_trans_factor=varying_trans_factor,
        if_wider_priors=if_wider_priors)

    if if_m_available:
        x_range = X_RANGE_WITH_M
        y_range = Y_RANGE_WITH_M
    else:
        x_range = X_RANGE_NO_M
        y_range = Y_RANGE_NO_M
    if varying_trans_factor:
        x_range = X_RANGE_VARYING_F
        y_range = Y_RANGE_VARYING_F

    plot_sa_for_varying_coverage(
        csv_file_name=csv_file_name,
        sim_duration=simulation_duration,
        fig_file_name=fig_file_name,
        x_range=x_range,
        y_range=y_range,
        interval=interval,
        varying_trans_factor=varying_trans_factor)


def plot_by_sens_spec(
        if_m_available, simulation_duration=None, calibration_seed=None):
    """
    plot cost-effectiveness plans varying sensitivity and specificity for this scenario for the availability of M
    :param if_m_available: (bool) if M is available for the first-line therapy
    :param simulation_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    """

    if if_m_available:
        x_range = X_RANGE_WITH_M
        y_range = Y_RANGE_WITH_M
    else:
        x_range = X_RANGE_NO_M
        y_range = Y_RANGE_NO_M

    for ab in ('CIP', 'TET'):

        for c in COVERAGE_VALUES:
            # plot a cost-effectiveness plan visualizing the performance of the rapid test for
            # the specified scenario for the availability of M, test coverage, and the antibiotic

            # get the filename of csv file where the scenario analysis and the
            # figure name which the figure for the scenario analysis should be saved as
            csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
                if_m_available_for_1st_tx=if_m_available,
                ab=ab,
                test_coverage=c,
                sim_duration=simulation_duration,
                calibration_seed=calibration_seed)

            plot_sa_for_specific_ab_and_coverage(
                csv_file_name=csv_file_name,
                fig_file_name=fig_file_name,
                ab=ab,
                include_sens_labels=True if ab == 'CIP' else False,
                test_coverage=c,
                x_range=x_range,
                y_range=y_range)

    for ab in ('CIP', 'TET'):

        #  plots a 3-panel figure of cost-effectiveness planes visualizing the performance of the rapid test for
        #  the specified scenario for the availability of M
        plot_by_sens_spec_coverage(
            ab=ab,
            include_sens_labels=True if ab == 'CIP' else False,
            if_m_available_for_1st_tx=if_m_available,
            sim_duration=simulation_duration,
            calibration_seed=calibration_seed
        )


if __name__ == "__main__":

    # plot for when sensitivity and specificity of tests have Beta distribution
    export_summary_and_plots_for_varying_coverage(
        if_m_available=True,  simulation_duration=SIM_DURATION, interval='c')

    export_summary_and_plots_for_varying_coverage(
        if_m_available=True, simulation_duration=35)

    export_summary_and_plots_for_varying_coverage(
        if_m_available=True, simulation_duration=SIM_DURATION, varying_trans_factor=True)

    # with wider prior distributions
    # for m in (True, False):
    #     export_summary_and_plots_for_varying_coverage(
    #         if_m_available=m, simulation_duration=SIM_DURATION, if_wider_priors=True)

    export_summary_and_plots_for_varying_coverage(
        if_m_available=True, simulation_duration=SIM_DURATION, calibration_seed=1)

    export_summary_and_plots_for_varying_coverage(
        if_m_available=False, simulation_duration=SIM_DURATION)

    export_summary_and_plots_for_varying_coverage(
        if_m_available=False, simulation_duration=35)

    plot_by_sens_spec(
        if_m_available=True, simulation_duration=SIM_DURATION)
