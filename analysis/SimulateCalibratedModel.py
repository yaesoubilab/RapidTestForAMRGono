from model.ModelSettings import GonoSettings
from model.Support import simulate_calibrated_model

IF_M_AVAILABLE_FOR_FIRST_TX = True
N_OF_TRAJS_TO_USE_FOR_SIMULATION = 32   # number of trajectories with the highest likelihood to keep

if __name__ == "__main__":

    # get model settings
    sets = GonoSettings()
    sets.update_settings(sens=1, spec=1, prob_rapid_test=0.75, # base: (0, 1, 0)
                         if_m_available_for_1st_tx=IF_M_AVAILABLE_FOR_FIRST_TX)

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sample_seeds_by_weights=False,
                              if_run_in_parallel=True,
                              settings=sets)
