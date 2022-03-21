Prev = {
    'ALL': [],
    'MSM': [4.47, 3.58, 5.36],   # 2016
    'MSW': [],
    'MSMW': [],
}

Rate = {
    'ALL': [],
    'MSM': [6508.0],   # 2018
    'MSW': [],
    'MSMW': [],
}

PercSymp = {
    'ALL': [],
    'MSM': [67.9, 63.0, 72.9],
    'MSW': [],
    'MSMW': [],
}


def get_prev(group):
    first_row = [0, None, None, None]
    second_row = [1] + Prev[group]

    return [first_row, second_row]


def get_rate(group):
    first_row = [0, None, None, None]
    second_row = [1] + Rate[group]

    return [first_row, second_row]


def get_prec_symp(group):
    first_row = [0, None, None, None]
    second_row = [1] + PercSymp[group]

    return [first_row, second_row]


def get_resist_profile(group):

    # read csv file
    return None


Prevalence = [
    [0, None, None, None],
    [1, 4.47, 3.58, 5.36],   # 2016
]

GonorrheaRate = [
    [0, None],
    [1, 6508.0],   # 2018
]

PercSymptomatic = [
    [0, None, None, None],
    [1, 67.9, 63.0, 72.9],
]

PercResistProfile = {
    'CIP-NS, TET-S, CRO-S': [
        [0, None, None, None],
        [1, 0.4, 0.2, 0.7]
    ],
    'CIP-S, TET-NS, CRP-S': [
        [0, None, None, None],
        [1, 56.0, 44.8, 67.2]
    ],
    'CIP-S, TET-S, CRO-NS': [
        [0, None, None, None],
        [1, 0.0]
    ],
    'CIP-NS, TET-NS, CRO-S': [
        [0, None, None, None],
        [1, 15.8, 12.6, 18.9]
    ],
    'CIP-NS, TET-S, CRO-NS': [
        [0, None, None, None],
        [1, 0.0]
    ],
    'CIP-S, TET-NS, CRO-NS': [
        [0, None, None, None],
        [1, 0.0]
    ],
    'CIP-NS, TET-NS, CRO-NS': [
        [0, None, None, None],
        [1, 0.0]
    ],
    'CIP-S, TET-S, CRO-S': [
        [0, None, None, None],
        [1, 27.8, 22.2, 33.3]
    ],
}
