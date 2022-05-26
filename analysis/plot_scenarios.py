from definitions import COVERAGE_VALUES
from model.plots import plot_scenarios, plot_scenario_sa, \
    get_scenarios_csv_filename_and_fig_filename, export_summary_of_scenarios

IF_M_AVAILABLE_FOR_FIRST_TX = True
TEST_COVERAGE = 0.75
X_RANGE_WITH_M = [-0.1, 6.1]
Y_RANGE_WITH_M = [-500, 1500]
X_RANGE_NO_M = [-0.1, 6.1]
Y_RANGE_NO_M = [-4000, 500]


def plot_scenarios_for_a_test_coverage(if_m_available_for_1st_tx, test_coverage):

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available_for_1st_tx,
        test_coverage=test_coverage)

    if if_m_available_for_1st_tx:
        x_range = X_RANGE_WITH_M
        y_range = Y_RANGE_WITH_M
    else:
        x_range = X_RANGE_NO_M
        y_range = Y_RANGE_NO_M

    plot_scenarios(
        csv_file_name=csv_file_name,
        fig_file_name=fig_file_name,
        test_coverage=test_coverage,
        x_range=x_range,
        y_range=y_range)


def plot_scenarios_sa(if_m_available_for_1st_tx):

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available_for_1st_tx)

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


if __name__ == "__main__":

    # for scenarios where drug M is and is not available for 1st-line therapy
    for m in (True, False):

        # export summary
        export_summary_of_scenarios(if_m_available_for_1st_tx=m)

        # 3 panel figure
        plot_scenarios_sa(if_m_available_for_1st_tx=m)

        for c in COVERAGE_VALUES:
            # plot
            plot_scenarios_for_a_test_coverage(
                if_m_available_for_1st_tx=m,
                test_coverage=c)


