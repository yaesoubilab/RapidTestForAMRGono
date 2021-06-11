from scipy.stats import norm


Prevalence = [
    [1, 4.47, 3.58, 5.36]   # 2016
]

GonorrheaRate = [
    [1, 6508.0]   # 2018
]

PercSymptomatic = [
    [1, 67.9, 63.0, 72.9],
]


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
