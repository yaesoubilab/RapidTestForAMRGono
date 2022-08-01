from deampy.format_functions import format_interval
from deampy.random_variats import Beta

# [mean, stDev, min, max]
CIP_SENS_DIST = 0.98, 0.01, 0, 1
CIP_SPEC_DIST = 0.98, 0.01, 0, 1
TET_SENS_DIST = 0.80, 0.1, 0, 1
TET_SPEC_DIST = 0.95, 0.02, 0, 1


def print_intervals(name, mean_std_min_max):

    interval = Beta.get_uncertainty_interval(
        alpha=0.05,
        mean=mean_std_min_max[0],
        st_dev=mean_std_min_max[1],
        minimum=mean_std_min_max[2],
        maximum=mean_std_min_max[3])

    print(name, format_interval(interval=interval, format='%', deci=1))


print_intervals(name='CIP SENS', mean_std_min_max=CIP_SENS_DIST)
print_intervals(name='CIP SPEC', mean_std_min_max=CIP_SPEC_DIST)
print_intervals(name='TET SENS', mean_std_min_max=TET_SENS_DIST)
print_intervals(name='TET SPEC', mean_std_min_max=TET_SPEC_DIST)
