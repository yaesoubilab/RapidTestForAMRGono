from model.ModelSettings import GonoSettings
from model.Support import simulate_calibrated_model

IF_M_AVAILABLE_FOR_FIRST_TX = True
N_OF_TRAJS_TO_USE_FOR_SIMULATION = 160   # number of trajectories with the highest likelihood to keep
# sensitivity, specificity, coverage
SEN, SPE, COV = 0.75, 0.975, 0.75     # base: (0, 1, 0)

if __name__ == "__main__":

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=IF_M_AVAILABLE_FOR_FIRST_TX)
    sets.update_settings(sens=SEN, spec=SPE, prob_rapid_test=COV)

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sample_seeds_by_weights=False,
                              if_run_in_parallel=True,
                              settings=sets)
