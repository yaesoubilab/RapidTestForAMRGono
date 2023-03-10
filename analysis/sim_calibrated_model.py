import sys

from deampy.in_out_functions import make_directory
from definitions import SIM_DURATION, ROOT_DIR, TRANSMISSION_FACTOR_VALUES
from model.model_settings import GonoSettings
from model.support import simulate_calibrated_model

"""
To simulate different scenarios with respect to the characteristics of rapid tests 
(sensitivity, specificity, coverage) 
"""

N_OF_TRAJS = 2   # number of trajectories to simulate


def simulate_calibrated(dict_test_characts=None,
                        if_m_available=True, if_wider_priors=False, transmission_factor=1.0):
    """
    simulates trajectories from the calibrated model
    :param dict_test_characts: (dictionary) of test characteristics with keys
        cip_sens: (float) sensitivity of the rapid test for CIP susceptibility (default 0)
        cip_spec: (float) specificity of the rapid test for CIP susceptibility (default 1)
        tet_sens: (float) sensitivity of the rapid test for TET susceptibility (default 0)
        tet_spec: (float) specificity of the rapid test for TET susceptibility (default 1)
        coverage: (float) coverage of the rapid test (default 0)
    :param if_wider_priors: (bool) if wider priors should be used
    :param if_m_available: (bool) if drug M is available for 1st line therapy
    :param transmission_factor: (float) transmission factor
    """

    if dict_test_characts is None:
        test_text = 'No DST'
        dict_test_characts = {'cip_sens': 0, 'cip_spec': 1, 'tet_sens': 0, 'tet_spec': 1, 'coverage': 0}
    else:
        for p in ('cip_sens', 'cip_spec', 'tet_sens', 'tet_spec'):
            if p not in dict_test_characts:
                dict_test_characts[p] = None

        test_text = 'p({}, {}) q({}, {}) c{}'.format(
            dict_test_characts['cip_sens'],
            dict_test_characts['tet_sens'],
            dict_test_characts['cip_spec'],
            dict_test_characts['tet_spec'],
            dict_test_characts['coverage'])

    if transmission_factor == 1.0:
        trans_factor_text = ''
    else:
        trans_factor_text = ' f{}'.format(transmission_factor)

    if if_m_available:
        figure_filename = 'With M ' + test_text + trans_factor_text
    else:
        figure_filename = 'No M ' + test_text + trans_factor_text

    # get model settings
    sets = GonoSettings(if_m_available_for_1st_tx=if_m_available, if_wider_priors=if_wider_priors)
    sets.update_settings(cip_sens=dict_test_characts['cip_sens'], cip_spec=dict_test_characts['cip_spec'],
                         tet_sens=dict_test_characts['tet_sens'], tet_spec=dict_test_characts['tet_spec'],
                         prob_rapid_test=dict_test_characts['coverage'], transmission_factor=transmission_factor)

    print('\n --- '+figure_filename+' ---')
    simulate_calibrated_model(n_of_sims=N_OF_TRAJS,
                              sim_duration=SIM_DURATION,
                              sample_seeds_by_weights=False,
                              if_run_in_parallel=True,
                              settings=sets,
                              figure_filename=figure_filename)


if __name__ == "__main__":

    make_directory(ROOT_DIR + '/analysis/outputs/')
    sys.stdout = open(ROOT_DIR + '/analysis/outputs/summary_simulating_calibrated_models.txt', 'w')

    # base (M is available)
    simulate_calibrated(if_m_available=True)

    # M and rapid DST is available
    simulate_calibrated(if_m_available=True, dict_test_characts={'coverage': 0.75})

    # worst-case scenario (neither M nor rapid DST is available)
    simulate_calibrated(if_m_available=False)

    # M not available but rapid DST is available
    simulate_calibrated(if_m_available=False, dict_test_characts={'coverage': 0.75})

    # change in transmission factor
    simulate_calibrated(if_m_available=True, transmission_factor=TRANSMISSION_FACTOR_VALUES[0])
    simulate_calibrated(if_m_available=True, transmission_factor=TRANSMISSION_FACTOR_VALUES[-1])

    sys.stdout.close()
