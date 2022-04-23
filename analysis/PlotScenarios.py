from definitions import COVERAGE_VALUES
from model.Plots import plot_scenarios

IF_M_AVAILABLE_FOR_FIRST_TX = True
TEST_COVERAGE = 0.75
X_RANGE = [-0.1, 6.1]
Y_RANGE = [-2500, 2500]


def plot_scenarios_for_a_text_coverage(if_m_available_for_1st_tx, test_coverage):

    if if_m_available_for_1st_tx:
        fig_file_name = 'figures/SA-with M-coverage {:.2f}.png'.format(c)
        csv_file_name = 'outputs-with-M/scenarios/simulated_scenarios.csv'
    else:
        fig_file_name = 'figures/SA-no M-coverage {:.2f}.png'.format(c)
        csv_file_name = 'outputs-no-M/scenarios/simulated_scenarios.csv'

    plot_scenarios(
        csv_file_name=csv_file_name,
        fig_file_name=fig_file_name,
        test_coverage=test_coverage,
        x_range=X_RANGE,
        y_range=Y_RANGE)


if __name__ == "__main__":

    for m in (True, False):
        for c in COVERAGE_VALUES:
            plot_scenarios_for_a_text_coverage(
                if_m_available_for_1st_tx=m,
                test_coverage=c)


