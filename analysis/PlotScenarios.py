from definitions import COVERAGE_VALUES
from model.Plots import plot_scenarios, get_scenarios_csv_filename_and_fig_filename, export_summary_of_scenarios

IF_M_AVAILABLE_FOR_FIRST_TX = True
TEST_COVERAGE = 0.75
X_RANGE = [-0.1, 6.1]
Y_RANGE = [-2500, 2500]


def plot_scenarios_for_a_text_coverage(if_m_available_for_1st_tx, test_coverage):

    # get the filename of csv file where the scenario analysis and the
    # figure name which the figure for the scenario analysis should be saved as
    csv_file_name, fig_file_name = get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx=if_m_available_for_1st_tx,
        test_coverage=test_coverage)


    plot_scenarios(
        csv_file_name=csv_file_name,
        fig_file_name=fig_file_name,
        test_coverage=test_coverage,
        x_range=X_RANGE,
        y_range=Y_RANGE)


if __name__ == "__main__":

    for m in (True, False):

        # export summary
        export_summary_of_scenarios(if_m_available_for_1st_tx=m)

        for c in COVERAGE_VALUES:

            # plot
            plot_scenarios_for_a_text_coverage(
                if_m_available_for_1st_tx=m,
                test_coverage=c)


