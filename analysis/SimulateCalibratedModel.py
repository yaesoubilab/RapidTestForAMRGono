from model.ModelSettings import GonoSettings
from model.Support import simulate_calibrated_model

N_OF_TRAJS_TO_USE_FOR_SIMULATION = 160   # number of trajectories with the highest likelihood to keep


def simulate_calibrated(sen=0.0, spe=1.0, coverage=0.0, if_m_available=True): # base: (0, 1, 0)

    if if_m_available:
        figure_filename = 'Calibrated with M p{} q{} c{}.png'.format(
            sen, spe, coverage)
    else:
        figure_filename = 'Calibration no M p{} q{} c{}.png'.format(
            sen, spe, coverage)

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available)
    sets.update_settings(sens=sen, spec=spe, prob_rapid_test=coverage)

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS_TO_USE_FOR_SIMULATION,
                              sample_seeds_by_weights=False,
                              if_run_in_parallel=True,
                              settings=sets,
                              figure_filename=figure_filename)


if __name__ == "__main__":

    # base (M is available)
    simulate_calibrated(sen=0, spe=1, coverage=0, if_m_available=True)

    # M and rapid DST are available
    simulate_calibrated(sen=0.75, spe=0.975, coverage=0.75, if_m_available=True)

    # worse case (M is not available)
    simulate_calibrated(sen=0, spe=1, coverage=0, if_m_available=False)

    # worse case but rapid DST are available
    simulate_calibrated(sen=0.75, spe=0.975, coverage=0.75, if_m_available=False)
