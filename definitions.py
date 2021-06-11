import os
from enum import Enum

from scipy.stats import norm

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

SYMP_STATES = ['Symp', 'Asym']
SUSP_PROFILES = ['PEN-R', 'PEN/CFX-R', 'SUS']


class SymStat(Enum):
    SYMP = 0
    ASYM = 1


class AB(Enum):
    PEN = 0
    CFX = 1


class SuspProfile(Enum):
    # we put SUS last because its prevalence is calculated after calculating the prevalence
    # of all other drugs
    PEN = 0
    PEN_CFX = 1
    SUS = 2


class SympSuspProfiles:
    # to convert (symptom state, susceptibility profile) to an index and vice versa

    def __init__(self, n_symp_stats, n_susp_profiles):
        self.nSympStats = n_symp_stats
        self.nSuspProfiles = n_susp_profiles
        self.length = n_symp_stats * n_susp_profiles

    def get_row_index(self, symp_state, susp_profile):
        return self.nSuspProfiles * symp_state + susp_profile

    def get_age_group_and_profile(self, i):
        return int(i / self.nSuspProfiles), i % self.nSympStats

    def get_str_susp_profile(self, symp_state, susp_profile):
        return '{}-{}'.format(SYMP_STATES[symp_state], SUSP_PROFILES[susp_profile])


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