from model.model_settings import GonoSettings
from model.support import simulate_calibrated_model

N_OF_TRAJS = 2   # number of trajectories to simulate


def simulate_calibrated(sen=0.0, spe=1.0, coverage=0.0, if_m_available=True):
    """
    simulates trajectories from the calibrated model
    :param sen: (float) sensitivity of the rapid test (default 0)
    :param spe: (float) specificity of the rapid test (default 1)
    :param coverage: (float) coverage of the rapid test (default 0)
    :param if_m_available: (bool) if drug M is available for 1st line therapy
    """

    if if_m_available:
        figure_filename = 'Calibrated with M p{} q{} c{}.png'.format(
            sen, spe, coverage)
    else:
        figure_filename = 'Calibration no M p{} q{} c{}.png'.format(
            sen, spe, coverage)

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available)
    sets.update_settings(sens=sen, spec=spe, prob_rapid_test=coverage)

    simulate_calibrated_model(n_of_sims=N_OF_TRAJS,
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
