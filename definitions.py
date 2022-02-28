import os
from enum import Enum

import numpy as np
from scipy.stats import norm

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

SIM_DURATION = 25
END_OF_WARM_UP = 1 # 6
END_OF_CALIB = 6
N_BREAKS_SENSITIVITY = 11
N_BREAKS_SPECIFICITY = 1

SYMP_STATES = ['Symp', 'Asym']
# we put SUS last because its prevalence is calculated after
# calculating the prevalence of other drugs
REST_PROFILES = ['CIP', 'TET', 'CFX', 'CIP+TET', 'CIP+CFX', 'TET+CFX', 'CIP+TET+CFX', 'SUSP']
ANTIBIOTICS = ['CIP', 'TET', 'CFX']


class SympStat(Enum):
    SYMP = 0
    ASYM = 1


class AB(Enum):
    CIP = 0
    TET = 1
    CFX = 2


class RestProfile(Enum):
    # we put SUS last because its prevalence is calculated after calculating the prevalence
    # of other drugs
    CIP = 0
    TET = 1
    CFX = 2
    CIP_TET = 3
    CIP_CFX = 4
    TET_CFX = 5
    CIP_TET_CFX = 6
    SUS = 7


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
        elif p == RestProfile.CFX:
            next_p = RestProfile.CIP_CFX
        elif p == RestProfile.TET_CFX:
            next_p = RestProfile.CIP_TET_CFX
        else:
            reason_for_failure = TreatmentOutcome.INEFFECTIVE
    elif a == AB.TET:
        if p == RestProfile.SUS:
            next_p = RestProfile.TET
        elif p == RestProfile.CIP:
            next_p = RestProfile.CIP_TET
        elif p == RestProfile.CFX:
            next_p = RestProfile.TET_CFX
        elif p == RestProfile.CIP_CFX:
            next_p = RestProfile.CIP_TET_CFX
        else:
            reason_for_failure = TreatmentOutcome.INEFFECTIVE
    elif a == AB.CFX:
        if p == RestProfile.SUS:
            next_p = RestProfile.CFX
        elif p == RestProfile.CIP:
            next_p = RestProfile.CIP_CFX
        elif p == RestProfile.TET:
            next_p = RestProfile.TET_CFX
        elif p == RestProfile.CIP_TET:
            next_p = RestProfile.CIP_TET_CFX
        else:
            reason_for_failure = TreatmentOutcome.INEFFECTIVE
    else:
        raise Exception('Invalid value for antibiotic.')

    return next_p.value, reason_for_failure


class ConvertSympAndSuspAndAntiBio:
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
        return '{}-{}|Tx with {}'.format(SYMP_STATES[symp_state], REST_PROFILES[rest_profile], ANTIBIOTICS[antibiotic])


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


def get_list_sensitivity_specificity(n_breaks_sensitivity, n_breaks_specificity):

    values = []
    if n_breaks_sensitivity == 1:
        values_sen = [0]
    else:
        values_sen = np.linspace(0, 1, n_breaks_sensitivity)
    if n_breaks_specificity == 1:
        values_spe = [1]
    else:
        values_spe = np.linspace(0, 1, n_breaks_specificity)

    for sens in reversed(values_sen):
        for spec in values_spe:
            values.append([sens, spec])

    return values


def get_scenario_names(n_breaks_sensitivity, n_breaks_specificity):

    scenario_names = []
    values_p_q = get_list_sensitivity_specificity(n_breaks_sensitivity=n_breaks_sensitivity,
                                                  n_breaks_specificity=n_breaks_specificity)
    for v in values_p_q:
        scenario_names.append('(p={:.2f}, q={:.2f})'.format(v[0], v[1]))

    return scenario_names
