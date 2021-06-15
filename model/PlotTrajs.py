import apace.analysis.Trajectories as A
from definitions import RestProfile, SympStat, REST_PROFILES, ConvertSympAndSuspAndAntiBio
from model import Data as D


def plot(prev_multiplier=52, incd_multiplier=1,
         obs_prev_multiplier=1, obs_incd_multiplier=1,
         filename='trajectories.png'):
    """
    :param prev_multiplier: (int) to multiply the simulation time to convert it to year, week, or day.
    :param incd_multiplier: (int) to multiply the simulation period to covert it to year, week, or day.
    :param obs_prev_multiplier: (int) to multiply the prevalence survey time to convert it to year, week, or day.
    :param obs_incd_multiplier: (int) to multiply the incidence survey period to covert it to year, week, or day.
    :param filename: (string) filename to save the trajectories as
    :return:
    """

    sim_outcomes = A.SimOutcomeTrajectories(csv_directory='outputs/trajectories')

    # defaults
    A.TIME_0 = 0  # 2014
    A.X_RANGE = (0, 11)
    A.X_TICKS = [A.TIME_0, 5]  # x-axis ticks (min at 0 with interval of 5)
    A.X_LABEL = 'Year'  # x-axis label
    A.TRAJ_TRANSPARENCY = 0.25

    # plot information
    S = A.TrajPlotInfo(outcome_name='In: S',
                       title='Susceptible',
                       x_multiplier=prev_multiplier,
                       y_range=(0, 1500000))

    Is = []
    Fs = []
    covert_symp_susp = ConvertSympAndSuspAndAntiBio(
        n_symp_stats=len(SympStat), n_susp_profiles=len(RestProfile))
    i = 0
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            str_symp_susp = covert_symp_susp.get_str_symp_susp(symp_state=s, susp_profile=p)
            # Is
            Is.append(A.TrajPlotInfo(outcome_name='In: I ' + str_symp_susp,
                                     title='I ' + str_symp_susp,
                                     x_multiplier=prev_multiplier))
            # Fs: infectious compartments after treatment failure
            Fs.append(A.TrajPlotInfo(outcome_name='In: F ' + str_symp_susp,
                                     title='F ' + str_symp_susp,
                                     x_multiplier=prev_multiplier))
            # increment i
            i += 1

    validation_filename = 'figures/(validation) ' + filename

    list_plot_info = Is
    list_plot_info.extend(Fs)
    list_plot_info.extend([S])
    sim_outcomes.plot_multi_panel(n_rows=5, n_cols=4,
                                  list_plot_info=list_plot_info,
                                  figure_size=(7, 7),
                                  file_name=validation_filename)

    # ------------- Calibration Figure ---------------

    prev = A.TrajPlotInfo(outcome_name='Prevalence',
                          title='Prevalence (%)',
                          x_multiplier=obs_prev_multiplier, y_multiplier=100,
                          y_range=(0, 10),
                          calibration_info=A.CalibrationTargetPlotInfo(
                              rows_of_data=D.Prevalence)
                          )
    gono_rate = A.TrajPlotInfo(outcome_name='Rate of gonorrhea cases',
                               title='Rate of gonorrhea cases\n(Per 100,000 MSM population)',
                               x_multiplier=obs_incd_multiplier, y_multiplier=100000,
                               y_range=(0, 10000),
                               calibration_info=A.CalibrationTargetPlotInfo(
                                   rows_of_data=D.GonorrheaRate)
                               )

    perc_symp = A.TrajPlotInfo(outcome_name='Proportion of cases symptomatic',
                               title='Percent gonorrhea cases\nthat are symptomatic (%)',
                               x_multiplier=obs_incd_multiplier,
                               y_multiplier=100, y_range=(0, 100),
                               calibration_info=A.CalibrationTargetPlotInfo(
                                   rows_of_data=D.PercSymptomatic)
                               )

    perc_cases_by_rest_profile = []
    for p in range(len(REST_PROFILES) - 1):
        perc_cases_by_rest_profile.append(A.TrajPlotInfo(
            outcome_name='Proportion of cases resistant to '+REST_PROFILES[p],
            title='Proportion of cases \nresistant to {} (%)'.format(REST_PROFILES[p]),
            x_multiplier=obs_incd_multiplier,
            y_multiplier=100, y_range=(0, 100))
        )

    calibration_filename = 'figures/(calibration) ' + filename

    list_plot_info=[prev, gono_rate, perc_symp]
    list_plot_info.extend(perc_cases_by_rest_profile)
    sim_outcomes.plot_multi_panel(n_rows=2, n_cols=3,
                                  list_plot_info=list_plot_info,
                                  figure_size=(6, 4.5), show_subplot_labels=True,
                                  file_name=calibration_filename)
