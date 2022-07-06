import os
from enum import Enum

from scipy.stats import norm

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

SIM_DURATION = 25
END_OF_WARM_UP = 1
END_OF_CALIB = 15
RAPID_TEST_COVERAGE = 0.75

# for sensitivity analysis
SEN_VALUES = (0.5, 0.75, 1.0)    # sensitivity
SPE_VALUES = (0.95, 0.975, 1.0)  # specificity
COVERAGE_VALUES = (0.5, 0.75, 1.0)      # coverage


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


def get_list_sens_spec_coverage():
    """
    :return: (list) of [sensitivity, specificity, coverage of rapid test]
    """

    values = []
    for cov in COVERAGE_VALUES:
        for sens in reversed(SEN_VALUES):
            for spec in SPE_VALUES:
                values.append([sens, spec, cov])

    return values


def get_scenario_names():
    """
    :return: (list) of scenario names based on [sensitivity, specificity, coverage of rapid test]
    """

    scenario_names = ['Status quo (no rapid test)']
    values_p_q_c = get_list_sens_spec_coverage()
    for v in values_p_q_c:
        scenario_names.append('(p={:.3f}, q={:.3f}, c={:.3f})'.format(
            v[0], v[1], v[2]))

    return scenario_names


def get_scenario_name(if_m_available, sim_duration=None, calibration_seed=None):
    """
    :param if_m_available:
    :param sim_duration:
    :param calibration_seed:
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
        name += '-{}'.format(calibration_seed)

    return name
