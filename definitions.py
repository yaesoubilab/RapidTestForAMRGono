import os
from enum import Enum

from scipy.stats import norm

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

SYMP_STATES = ['Symp', 'Asym']
# we put SUS last because its prevalence is calculated after
# calculating the prevalence of all other drugs
REST_PROFILES = ['PEN-R', 'CFX-R', 'PEN+CFX-R', 'SUS']
ANTIBIOTICS = ['PEN', 'CFX']


class SympStat(Enum):
    SYMP = 0
    ASYM = 1


class AB(Enum):
    PEN = 0
    CFX = 1


class RestProfile(Enum):
    # we put SUS last because its prevalence is calculated after calculating the prevalence
    # of all other drugs
    PEN = 0
    CFX = 1
    PEN_CFX = 2
    SUS = 3


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

    def get_str_symp_susp(self, symp_state, susp_profile):
        return '{}-{}'.format(SYMP_STATES[symp_state], REST_PROFILES[susp_profile])

    def get_str_symp_rest_antibio(self, symp_state, rest_profile, antibiotic):
        return '{}-{}|Tx with {}'.format(SYMP_STATES[symp_state], REST_PROFILES[rest_profile], ANTIBIOTICS[antibiotic])


def get_survey_size(mean, l, u, multiplier=1):

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