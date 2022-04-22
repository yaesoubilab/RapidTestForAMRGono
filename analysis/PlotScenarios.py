from definitions import get_scenario_names
from model.Plots import plot_scenarios


IF_M_AVAILABLE_FOR_FIRST_TX = False
X_RANGE = [-0.1, 8.1]
Y_RANGE = [-6000, 5000]


if __name__ == "__main__":

    if IF_M_AVAILABLE_FOR_FIRST_TX:
        fig_file_name = 'figures/Changing specificity-with-M.png'
        csv_file_name = 'outputs-with-M/scenarios/simulated_scenarios.csv'
    else:
        fig_file_name = 'figures/Changing specificity-no-M.png'
        csv_file_name = 'outputs-no-M/scenarios/simulated_scenarios.csv'

    plot_scenarios(
        scenario_names=get_scenario_names(),
        csv_file_name=csv_file_name,
        fig_file_name=fig_file_name,
        x_range=X_RANGE,
        y_range=Y_RANGE)


