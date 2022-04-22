from model.Plots import plot_scenarios


IF_M_AVAILABLE_FOR_FIRST_TX = True
TEST_COVERAGE = 0.75
X_RANGE = [-0.1, 8.1]
Y_RANGE = [-6000, 5000]


if __name__ == "__main__":

    if IF_M_AVAILABLE_FOR_FIRST_TX:
        fig_file_name = 'figures/SA-with M-coverage {}.png'.format(TEST_COVERAGE)
        csv_file_name = 'outputs-with-M/scenarios/simulated_scenarios.csv'
    else:
        fig_file_name = 'figures/SA-no M-coverage {}.png'.format(TEST_COVERAGE)
        csv_file_name = 'outputs-no-M/scenarios/simulated_scenarios.csv'

    plot_scenarios(
        csv_file_name=csv_file_name,
        fig_file_name=fig_file_name,
        test_coverage=TEST_COVERAGE,
        x_range=X_RANGE,
        y_range=Y_RANGE)


