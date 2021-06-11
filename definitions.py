import os
from enum import Enum

from scipy.stats import norm

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

SUSP_PROFILES = ['SUS', 'PEN-R', 'PEN/CFX-R']


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


class AgeGroupsProfiles:
    # to convert (age group index, profile index) to an index and vice versa

    def __init__(self, n_age_groups, n_profiles):
        self.nAgeGroups = n_age_groups
        self.nProfiles = n_profiles
        self.length = n_age_groups * n_profiles

    def get_row_index(self, age_group, profile):
        return self.nProfiles * age_group + profile

    def get_age_group_and_profile(self, i):
        return int(i/self.nProfiles), i % self.nAgeGroups

    def get_str_age_profile(self, age_group, profile):

        if profile == 0:
            return AGES[age_group] + '-Current'
        elif profile == 1:
            return AGES[age_group] + '-Novel'
        else:
            raise ValueError
        #return 'age/profile ({},{})'.format(age_group, profile)

    def get_str_age(self, age_group):
        return AGES[age_group]
        # return 'age {}'.format(age_group)


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