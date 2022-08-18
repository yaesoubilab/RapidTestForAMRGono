import sys

from deampy.in_out_functions import make_directory
from definitions import SIM_DURATION, ROOT_DIR
from model.model_settings import GonoSettings
from model.support import simulate_calibrated_model

"""
To simulate different scenarios with respect to the characteristics of rapid tests 
(sensitivity, specificity, coverage) 
"""

N_OF_TRAJS = 32   # number of trajectories to simulate


def simulate_calibrated(cip_sens=None, cip_spec=None,
                        tet_sens=None, tet_spec=None, coverage=0.0, if_m_available=True):
    """
    simulates trajectories from the calibrated model
    :param cip_sens: (float) sensitivity of the rapid test for CIP susceptibility (default 0)
    :param cip_spec: (float) specificity of the rapid test for CIP susceptibility (default 1)
    :param tet_sens: (float) sensitivity of the rapid test for TET susceptibility (default 0)
    :param tet_spec: (float) specificity of the rapid test for TET susceptibility (default 1)
    :param coverage: (float) coverage of the rapid test (default 0)
    :param if_m_available: (bool) if drug M is available for 1st line therapy
    """

    if if_m_available:
        figure_filename = 'With M p({}, {}) q({}, {}) c{}.png'.format(
            cip_sens, tet_sens, cip_spec, tet_spec, coverage)
    else:
        figure_filename = 'No M p({}, {}) q({}, {}) c{}.png'.format(
            cip_sens, tet_sens, cip_spec, tet_spec, coverage)

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available)
    sets.update_settings(cip_sens=cip_sens, cip_spec=cip_spec,
                         tet_sens=tet_sens, tet_spec=tet_spec,
                         prob_rapid_test=coverage)

    print('\n --- '+figure_filename+' ---')
    simulate_calibrated_model(n_of_sims=N_OF_TRAJS,
                              sim_duration=SIM_DURATION,
                              sample_seeds_by_weights=False,
                              if_run_in_parallel=True,
                              settings=sets,
                              figure_filename=figure_filename)


if __name__ == "__main__":

    make_directory(ROOT_DIR + '/analysis/outputs/')
    sys.stdout = open(ROOT_DIR + '/analysis/outputs/scenario_performance.txt', 'w')

    # base (M is available)
    simulate_calibrated(cip_sens=0, cip_spec=1,
                        tet_sens=0, tet_spec=1,
                        coverage=0, if_m_available=True)
    #
    # M and rapid DST is available
    simulate_calibrated(cip_sens=None, cip_spec=None,
                        tet_sens=None, tet_spec=None,
                        coverage=0.75, if_m_available=True)

    # worst-case scenario (neither M nor rapid DST is available)
    simulate_calibrated(cip_sens=0, cip_spec=1,
                        tet_sens=0, tet_spec=1,
                        coverage=0, if_m_available=False)

    # M not available but and rapid DST is available
    simulate_calibrated(cip_sens=None, cip_spec=None,
                        tet_sens=None, tet_spec=None,
                        coverage=0.75, if_m_available=False)

    sys.stdout.close()
