from model.support import simulate_multi_trajectories

if __name__ == "__main__":

    simulate_multi_trajectories(n=32, sim_duration=25, figure_filename='Uncalibrated.png', if_run_in_parallel=True)
