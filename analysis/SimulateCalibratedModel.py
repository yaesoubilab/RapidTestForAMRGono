from model.ModelSettings import GonoSettings
from model.Support import simulate_calibrated_model

N_OF_TRAJS_TO_USE_FOR_SIMULATION = 64   # number of trajectories with the highest likelihood to keep

if __name__ == "__main__":

    # get model settings
    sets = GonoSettings()
    sets.update_settings(sens=0, spec=1, prob_rapid_test=0)  # base: (0, 1, 0)

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sample_seeds_by_weights=False,
                              if_run_in_parallel=True,
                              settings=sets)
