from model.Support import simulate_calibrated_model

N_OF_TRAJS_TO_USE_FOR_SIMULATION = 50   # number of trajectories with the highest likelihood to keep
IF_DISRUPTION = False

if __name__ == "__main__":

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sample_seeds_by_weights=False)
