import os
from enum import Enum

from scipy.stats import norm

FIG_EXT = 'pdf'  # 'png'
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

SIM_DURATION = 25
END_OF_WARM_UP = 1
END_OF_CALIB = 15
RAPID_TEST_COVERAGE = 0.75

# for sensitivity analysis
# mean, st_dev, min, max of beta distributions
CIP_SENS_DIST = 0.98, 0.01, 0, 1
CIP_SPEC_DIST = 0.98, 0.01, 0, 1
TET_SENS_DIST = 0.855, 0.014, 0, 1
TET_SPEC_DIST = 0.965, 0.008, 0, 1

CIP_SENS_VALUES = (0.5, 0.75, 1.0)
CIP_SPEC_VALUES = (0.9, 0.95, 1.0)
TET_SENS_VALUES = (0.5, 0.75, 1.0)
TET_SPEC_VALUES = (0.9, 0.95, 1.0)

COVERAGE_VALUES = (0.5, 0.75, 1.0)      # coverage
TRANSMISSION_FACTOR_VALUES = (0.98, 1.0, 1.02)
TRANSMISSION_FACTOR_COLORS = ('red', 'green', 'blue')

# plots
COLOR_VARYING_COVERAGE = 'blue'
COLOR_BY_SPEC = ['red', 'green', 'purple']
SINGLE_FIG_SIZE = (4.6, 4.6)
TRIPLE_FIG_SIZE = (8, 3.5)

EFFECT_OUTCOME = 'Time-averaged proportion of cases treated successfully with CIP, TET, or CRO ' \
                 '(average incidence after epidemic warm-up)'
COST_OUTCOME = 'Rate of gonorrhea cases (average incidence after epidemic warm-up)'
EFFECT_COST_LABELS = ('Change in the combined effective lifespan of\n first-line antibiotics (years)',
                      'Change in the annual rate of gonorrhea\n(per 100,000 MSM population)')
EFFECT_COST_LABELS_NO_LINE_BREAK = ('Change in the combined effective lifespan of first-line antibiotics (years)',
                                    'Change in the annual rate of gonorrhea\n(per 100,000 MSM population)')


SYMP_STATES = ['Symp', 'Asym']
REST_PROFILES = ['CIP-S, TET-S, CRO-S',        # SUSP
                 'CIP-NS, TET-S, CRO-S',       # not susceptible to CIP
                 'CIP-S, TET-NS, CRO-S',       # not susceptible to TET
                 'CIP-S, TET-S, CRO-NS',       # not susceptible to CRO
                 'CIP-NS, TET-NS, CRO-S',      # not susceptible to CIP+TET
                 'CIP-NS, TET-S, CRO-NS',      # not susceptible to CIP+CRO
                 'CIP-S, TET-NS, CRO-NS',      # not susceptible to TET+CRO
                 'CIP-NS, TET-NS, CRO-NS',     # not susceptible to CIP+TET+CRO
                 ]

ANTIBIOTICS = ['CIP', 'TET', 'CRO']


class SympStat(Enum):
    SYMP = 0
    ASYM = 1


class AB(Enum):
    CIP = 0
    TET = 1
    CRO = 2


class RestProfile(Enum):
    SUS = 0
    CIP = 1
    TET = 2
    CRO = 3
    CIP_TET = 4
    CIP_CRO = 5
    TET_CRO = 6
    CIP_TET_CRO = 7


class TreatmentOutcome(Enum):
    SUCCESS = 0
    RESISTANCE = 1
    INEFFECTIVE = 2


def get_profile_after_resit_or_failure(rest_profile, antibiotic):
    """
    :param rest_profile: (int) index of current resistance profile
    :param antibiotic: (int) index of antibiotic to use
    :return: (resistance profile after treatment, reason for treatment failure)
    """

    p = RestProfile(rest_profile)
    a = AB(antibiotic)
    # default values
    next_p = p
    reason_for_failure = TreatmentOutcome.RESISTANCE

    if a == AB.CIP:
        if p == RestProfile.SUS:
            next_p = RestProfile.CIP
        elif p == RestProfile.TET:
            next_p = RestProfile.CIP_TET
        elif p == RestProfile.CRO:
            next_p = RestProfile.CIP_CRO
        elif p == RestProfile.TET_CRO:
            next_p = RestProfile.CIP_TET_CRO
        else:
            reason_for_failure = TreatmentOutcome.INEFFECTIVE
    elif a == AB.TET:
        if p == RestProfile.SUS:
            next_p = RestProfile.TET
        elif p == RestProfile.CIP:
            next_p = RestProfile.CIP_TET
        elif p == RestProfile.CRO:
            next_p = RestProfile.TET_CRO
        elif p == RestProfile.CIP_CRO:
            next_p = RestProfile.CIP_TET_CRO
        else:
            reason_for_failure = TreatmentOutcome.INEFFECTIVE
    elif a == AB.CRO:
        if p == RestProfile.SUS:
            next_p = RestProfile.CRO
        elif p == RestProfile.CIP:
            next_p = RestProfile.CIP_CRO
        elif p == RestProfile.TET:
            next_p = RestProfile.TET_CRO
        elif p == RestProfile.CIP_TET:
            next_p = RestProfile.CIP_TET_CRO
        else:
            reason_for_failure = TreatmentOutcome.INEFFECTIVE
    else:
        raise Exception('Invalid value for antibiotic.')

    return next_p.value, reason_for_failure


class ConvertSympAndResitAndAntiBio:
    # to convert (symptom state, resistance profile, antibiotic) to an index and vice versa

    def __init__(self, n_symp_stats, n_rest_profiles, n_antibiotics=None):
        self.nSympStats = n_symp_stats
        self.nRestProfiles = n_rest_profiles
        self.nAntiBiotics = n_antibiotics
        if n_antibiotics in (None, 0):
            self.length = n_symp_stats * n_rest_profiles
        else:
            self.length = n_symp_stats * n_rest_profiles * n_antibiotics

    def get_row_index(self, symp_state, rest_profile, antibiotic=None):
        if self.nAntiBiotics in (None, 0):
            return self.nRestProfiles * symp_state + rest_profile
        else:
            return (self.nRestProfiles * self.nAntiBiotics) * symp_state \
                   + self.nAntiBiotics * rest_profile + antibiotic

    def get_symp_and_profile(self, i):
        if self.nAntiBiotics in (None, 0):
            return int(i / self.nRestProfiles), i % self.nSympStats
        else:
            return None

    def get_str_susp(self, susp_profile):
        return REST_PROFILES[susp_profile]

    def get_str_symp_susp(self, symp_state, rest_profile):
        return '{}-{}'.format(SYMP_STATES[symp_state], REST_PROFILES[rest_profile])

    def get_str_symp_rest_antibio(self, symp_state, rest_profile, antibiotic):
        return 'Tx-{} in {}-{}'.format(
            ANTIBIOTICS[antibiotic], SYMP_STATES[symp_state], REST_PROFILES[rest_profile])


def get_survey_size(mean, l, u, multiplier=1):

    if mean is None:
        return None

    mean *= multiplier
    l *= multiplier
    u *= multiplier

    # st-dev
    var = mean*(1-mean)

    # half-width
    hw = (u-l)/2

    # z
    z = norm.ppf(1-0.05/2)

    return var * pow(z/hw, 2)


def get_name_of_scenario_analysis(cip_sens, cip_spec, tet_sens, tet_spec, coverage, transmission_factor):
    """
    :return: the name of this scenario
    """

    cip_sens_name = 'None' if cip_sens is None else '{:.3f}'.format(cip_sens)
    cip_spec_name = 'None' if cip_spec is None else '{:.3f}'.format(cip_spec)
    tet_sens_name = 'None' if tet_sens is None else '{:.3f}'.format(tet_sens)
    tet_spec_name = 'None' if tet_spec is None else '{:.3f}'.format(tet_spec)

    return 'p=({}, {}), q=({}, {}), c={:0.3f}, f={:0.2f})'.format(
        cip_sens_name, tet_sens_name, cip_spec_name, tet_spec_name, coverage, transmission_factor)


def get_sens_analysis_names_and_definitions(vary_sens_spec=False, vary_transm_factor=False):
    """
    :param vary_sens_spec: (bool) if include analyses to vary sensitivity and specificity
        of CIP and TET
    :param vary_transm_factor: (bool) if includes analyses to vary transmission factor
    :return: (tuple of two list) (scenario names, scenarios definitions)
        according to:
        [CIP-sensitivity, CIP-specificity,
        TET-sensitivity, TET-specificity,
        coverage of rapid test, transmission factor]
    """

    names = ['Status quo (no rapid test)']
    definitions = [[0.0, 1.0, 0.0, 1.0, 0.0, 1.0]]

    if not vary_transm_factor:
        # varying coverage but using beta distributions for sensitivity and specificity of tests
        for cov in COVERAGE_VALUES:
            names.append(get_name_of_scenario_analysis(
                cip_sens=None, cip_spec=None, tet_sens=None, tet_spec=None,
                coverage=cov, transmission_factor=1.0))
            definitions.append([None, None, None, None, cov, 1.0])

    else:
        # varying transmission factor and coverage but
        # using beta distributions for sensitivity and specificity of tests
        for f in TRANSMISSION_FACTOR_VALUES:
            for cov in COVERAGE_VALUES:
                names.append(get_name_of_scenario_analysis(
                    cip_sens=None, cip_spec=None, tet_sens=None, tet_spec=None,
                    coverage=cov, transmission_factor=f))
                definitions.append([None, None, None, None, cov, f])

    if vary_sens_spec:
        # use beta distributions for CIP and vary characteristics of TET
        for cov in COVERAGE_VALUES:
            for tet_sens in reversed(TET_SENS_VALUES):
                for tet_spec in TET_SPEC_VALUES:
                    names.append(get_name_of_scenario_analysis(
                        cip_sens=None, cip_spec=None, tet_sens=tet_sens, tet_spec=tet_spec,
                        coverage=cov, transmission_factor=1.0))
                    definitions.append([None, None, tet_sens, tet_spec, cov, 1.0])

        # use beta distributions for TET and vary characteristics of CIP
        for cov in COVERAGE_VALUES:
            for cip_sens in reversed(CIP_SENS_VALUES):
                for cip_spec in CIP_SPEC_VALUES:
                    names.append(get_name_of_scenario_analysis(
                        cip_sens=cip_sens, cip_spec=cip_spec, tet_sens=None, tet_spec=None,
                        coverage=cov, transmission_factor=1))
                    definitions.append([cip_sens, cip_spec, None, None, cov, 1])

    return names, definitions


def get_scenario_name(if_m_available, sim_duration=None, calibration_seed=None,
                      if_wider_priors=False, if_varying_transmission_factor=False):
    """
    :param if_m_available: (bool) if M is available for first-line therapy
    :param sim_duration: (float) simulation duration
    :param calibration_seed: (int or None) calibration seed (for sensitivity analysis)
    :param if_wider_priors: (bool) set to True if the winder prior distributions should be used
    :param if_varying_transmission_factor: (bool) set to True if transmission factor will be varied
        (for sensitivity analysis)
    :return: the name the scenario being simulated
    """

    name = ''

    if if_m_available:
        name += 'with M'
    else:
        name += 'no M'

    if sim_duration is not None:
        name += '-{}yrs'.format(sim_duration)

    if calibration_seed is not None:
        name += '-seed{}'.format(calibration_seed)

    if if_wider_priors:
        name += '-wider priors'

    if if_varying_transmission_factor:
        name += '-varying f'

    return name


def get_traj_fig_name(if_m_available, dict_test_characts, transmission_factor):

    if dict_test_characts is None:
        test_text = 'No DST'
    else:
        test_p_q_provided = True
        for p in ('cip_sens', 'cip_spec', 'tet_sens', 'tet_spec'):
            if p not in dict_test_characts:
                dict_test_characts[p] = None
                test_p_q_provided = False

        if test_p_q_provided:
            test_text = 'p({}, {}) q({}, {}) c{}'.format(
                dict_test_characts['cip_sens'],
                dict_test_characts['tet_sens'],
                dict_test_characts['cip_spec'],
                dict_test_characts['tet_spec'],
                int(dict_test_characts['coverage']*100))
        else:
            test_text = 'With DST c{}'.format(
                int(dict_test_characts['coverage'] * 100))

    if transmission_factor == 1.0:
        trans_factor_text = ''
    else:
        trans_factor_text = ' f{}'.format(int(100*transmission_factor))

    if if_m_available:
        figure_filename = 'With M ' + test_text + trans_factor_text
    else:
        figure_filename = 'No M ' + test_text + trans_factor_text

    return '{}.{}'.format(figure_filename, FIG_EXT)
