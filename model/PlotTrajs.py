import apace.analysis.Trajectories as A
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
    A.TRAJ_TRANSPARENCY = 0.2

    # plot information
    S = A.TrajPlotInfo(outcome_name='In: S',
                       title='Susceptible',
                       x_multiplier=prev_multiplier,
                       y_range=(0, 1500000))
    I0_Sym = A.TrajPlotInfo(outcome_name='In: I | Sus | Sym',
                            title='I | Sus | Sym',
                            x_multiplier=prev_multiplier)
    I0_Asym = A.TrajPlotInfo(outcome_name='In: I | Sus | Asym',
                             title='I | Sus | Asym',
                             x_multiplier=prev_multiplier)
    IA_Sym = A.TrajPlotInfo(outcome_name='In: I | Res | Sym',
                            title='I | Res | Sym',
                            x_multiplier=prev_multiplier)
    IA_Asym = A.TrajPlotInfo(outcome_name='In: I | Res | Asym',
                             title='I | Res | Asym',
                             x_multiplier=prev_multiplier)
    W2_IA_Sym = A.TrajPlotInfo(outcome_name='In: Waiting 2nd | I | Res | Sym',
                               title='Waiting 2nd | I | Res | Sym',
                               x_multiplier=prev_multiplier)

    To_I0_Sym = A.TrajPlotInfo(outcome_name='To: I | Sus | Sym',
                               title='To: I | Sus | Sym',
                               x_multiplier=incd_multiplier)
    To_I0_Asym = A.TrajPlotInfo(outcome_name='To: I | Sus | Asym',
                                title='To: I | Sus | Asym',
                                x_multiplier=incd_multiplier)
    To_IA_Sym = A.TrajPlotInfo(outcome_name='To: I | Res | Sym',
                               title='To: I | Res | Sym',
                               x_multiplier=incd_multiplier)
    To_IA_Asym = A.TrajPlotInfo(outcome_name='To: I | Res | Asym',
                                title='To: I | Res | Asym',
                                x_multiplier=incd_multiplier)

    validation_filename = 'figures/(validation) ' + filename

    sim_outcomes.plot_multi_panel(n_rows=2, n_cols=5,
                                  list_plot_info=[S, I0_Sym, I0_Asym, IA_Sym, IA_Asym,
                                                  W2_IA_Sym, To_I0_Sym, To_I0_Asym, To_IA_Sym, To_IA_Asym],
                                  figure_size=(9, 4),
                                  file_name=validation_filename)

    # ------------- Calibration Figure ---------------

    Prev = A.TrajPlotInfo(outcome_name='Prevalence',
                          title='Prevalence (%)',
                          x_multiplier=obs_prev_multiplier, y_multiplier=100,
                          y_range=(0, 10),
                          calibration_info=A.CalibrationTargetPlotInfo(
                              rows_of_data=D.Prevalence)
                          )
    GonoRate = A.TrajPlotInfo(outcome_name='Rate of gonorrhea cases',
                              title='Rate of gonorrhea cases\n(Per 100,000 MSM population)',
                              x_multiplier=obs_incd_multiplier, y_multiplier=100000,
                              y_range=(0, 10000),
                              calibration_info=A.CalibrationTargetPlotInfo(
                                  rows_of_data=D.GonorrheaRate)
                              )

    PercentSymptomatic = A.TrajPlotInfo(outcome_name='Proportion of cases symptomatic',
                                        title='Percent gonorrhea cases\nthat are symptomatic (%)',
                                        x_multiplier=obs_incd_multiplier,
                                        y_multiplier=100, y_range=(0, 100),
                                        calibration_info=A.CalibrationTargetPlotInfo(
                                            rows_of_data=D.PercSymptomatic)
                                        )

    PercentCasesResistantA = A.TrajPlotInfo(outcome_name='Proportion of cases resistant',
                                            title='Proportion of cases with\nresistant strain (%)',
                                            x_multiplier=obs_incd_multiplier,
                                            y_multiplier=100, y_range=(0, 100))

    calibration_filename = 'figures/(calibration) ' + filename

    sim_outcomes.plot_multi_panel(n_rows=2, n_cols=2,
                                  list_plot_info=[Prev, GonoRate, PercentSymptomatic, PercentCasesResistantA],
                                  figure_size=(4.5, 4.5), show_subplot_labels=True,
                                  file_name=calibration_filename)
