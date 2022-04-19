import numpy as np

import apace.analysis.Scenarios as S
import apace.analysis.Trajectories as A
import apace.analysis.VisualizeScenarios as V
from definitions import RestProfile, SympStat, REST_PROFILES, \
    ConvertSympAndResitAndAntiBio, SIM_DURATION, ANTIBIOTICS, N_BREAKS_SPECIFICITY, MIN_SPE
from model import Data as D

A.SUBPLOT_W_SPACE = 0.25
S.POLY_DEGREES = 1
COLORS = ['purple', 'blue', 'red', 'green', 'orange', 'brown']


def plot_trajectories(prev_multiplier=52, incd_multiplier=1,
                      obs_prev_multiplier=1, obs_incd_multiplier=1,
                      dir_of_trajs='outputs/trajectories',
                      filename='trajectories.png'):
    """
    :param prev_multiplier: (int) to multiply the simulation time to convert it to year, week, or day.
    :param incd_multiplier: (int) to multiply the simulation period to covert it to year, week, or day.
    :param obs_prev_multiplier: (int) to multiply the prevalence survey time to convert it to year, week, or day.
    :param obs_incd_multiplier: (int) to multiply the incidence survey period to covert it to year, week, or day.
    :param dir_of_trajs: (string) directory where simulated trajectories are located
    :param filename: (string) filename to save the trajectories as
    :return:
    """

    sim_outcomes = A.SimOutcomeTrajectories(csv_directory=dir_of_trajs)

    # defaults
    A.TIME_0 = 0  # 2014
    A.X_RANGE = (0, SIM_DURATION+1)
    A.X_TICKS = [A.TIME_0, 5]  # x-axis ticks (min at 0 with interval of 5)
    A.X_LABEL = 'Year'  # x-axis label
    A.TRAJ_TRANSPARENCY = 0.75 if filename == 'onetraj' else 0.25

    # plot information
    S = A.TrajPlotInfo(outcome_name='In: S',
                       title='Susceptible',
                       x_multiplier=prev_multiplier,
                       y_range=(0, 1500000))
    pop = A.TrajPlotInfo(outcome_name='Population size',
                         title='Population',
                         x_multiplier=prev_multiplier,
                         y_range=(0.95*pow(10, 6), 1.05*pow(10, 6)))

    Is = []
    Fs = []
    covert_symp_susp = ConvertSympAndResitAndAntiBio(
        n_symp_stats=len(SympStat), n_rest_profiles=len(RestProfile))
    i = 0
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            str_symp_susp = covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p)
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

    sim_outcomes.plot_multi_panel(n_rows=4, n_cols=4,
                                  list_plot_info=Is,
                                  figure_size=(7, 7),
                                  file_name='figures/(valid-Is) ' + filename)
    sim_outcomes.plot_multi_panel(n_rows=4, n_cols=4,
                                  list_plot_info=Fs,
                                  figure_size=(7, 7),
                                  file_name='figures/(valid-Fs) ' + filename)

    # sim_outcomes.plot_multi_panel(n_rows=1, n_cols=2,
    #                               list_plot_info=[S, pop],
    #                               figure_size=(2*2.2, 1*2.2),
    #                               file_name='figures/(valid-S) ' + filename)

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
                               y_range=(0, 12500),
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
    for p in range(len(REST_PROFILES)):
        perc_cases_by_rest_profile.append(A.TrajPlotInfo(
            outcome_name='Proportion of cases resistant to '+REST_PROFILES[p],
            title='Cases with profile\n{} (%)'.format(REST_PROFILES[p]),
            x_multiplier=obs_incd_multiplier,
            y_multiplier=100, y_range=(0, 100),
            calibration_info=A.CalibrationTargetPlotInfo(
                rows_of_data=D.PercResistProfile[REST_PROFILES[p]]))
        )

    calibration_filename = 'figures/(summary) ' + filename

    list_plot_info = [prev, gono_rate, perc_symp]
    list_plot_info.extend(perc_cases_by_rest_profile)
    sim_outcomes.plot_multi_panel(n_rows=4, n_cols=3,
                                  list_plot_info=list_plot_info,
                                  figure_size=(3*2.2, 4*2.3), show_subplot_labels=True,
                                  file_name=calibration_filename)

    # ------------- Successful treatment with different antibiotics ---------------
    Txs = []
    for a in ANTIBIOTICS:
        Txs.append(A.TrajPlotInfo(outcome_name='To: Counting success tx with ' + a,
                                  title='Successful Tx-' + a,
                                  x_multiplier=incd_multiplier))
    Txs.append(A.TrajPlotInfo(outcome_name='To: 1st-Tx with M',
                              title='Successful 1st-Tx-M',
                              x_multiplier=incd_multiplier))
    Txs.append(A.TrajPlotInfo(outcome_name='To: Tx with M',
                              title='Successful Tx-M',
                              x_multiplier=incd_multiplier))
    sim_outcomes.plot_multi_panel(n_rows=2, n_cols=3,
                                  list_plot_info=Txs,
                                  figure_size=(3*2.2, 2*2.2), show_subplot_labels=True,
                                  file_name='figures/(valid-Tx) ' + filename)


def plot_scenarios(scenario_names, csv_file_name, fig_file_name, x_range, y_range):

    # read scenarios into a dataframe
    scenarios_df = S.ScenarioDataFrame(csv_file_name=csv_file_name)

    # get a specific outcome from a specific scenario
    print('\nScenario name: ', 'Rate of gonorrhea cases | Effective lifespan of CIP, TET, and CFX')
    for name in scenario_names:
        rate = scenarios_df.get_mean_interval(
            scenario_name=name,
            outcome_name='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
            interval_type='c', deci=3)
        life = scenarios_df.get_mean_interval(
            scenario_name=name,
            outcome_name='Proportion of cases treated with CIP, TET, or CRO (average incidence after epidemic warm-up)',
            interval_type='c', deci=3, multiplier=SIM_DURATION)
        print('{}: \t{} | {}'.format(name, rate, life))

    # plot CEA
    S.ERROR_BAR_ALPHA = 0.2
    # # read scenarios into a dataframe
    # df_scenarios = S.ScenarioDataFrame(
    #     csv_file_name='outputs/scenarios/simulated_scenarios.csv')

    # sets of scenarios to display on the cost-effectiveness plain
    list_of_scenario_sets = []
    list_if_remove_base_scenario = [True] * N_BREAKS_SPECIFICITY
    for i, spec in enumerate(np.linspace(MIN_SPE, 1, N_BREAKS_SPECIFICITY)):
        list_of_scenario_sets.append(S.SetOfScenarios(
            name='Specificity = {:.2f}'.format(spec),
            scenario_df=scenarios_df,
            color=COLORS[i],
            marker='o',
            conditions_on_variables=[
                S.ConditionOnVariable(var_name='sensitivity', if_included_in_label=True, label_format='{:.1f}'),
                S.ConditionOnVariable(var_name='specificity', values=[spec])],
            if_find_frontier=False,
            if_show_fitted_curve=True,
            labels_shift_x=0.03,
            labels_shift_y=0.00)
        )

    V.plot_sets_of_scenarios(
        list_of_scenario_sets=list_of_scenario_sets,
        name_of_base_scenario='Status quo (no rapid test)',
        list_if_remove_base_scenario=list_if_remove_base_scenario,
        effect_outcome='Proportion of cases treated with CIP, TET, or CRO (average incidence after epidemic warm-up)',
        cost_outcome='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
        labels=('Change in the effective lifespan of\nCIP, TET, and CRO',
                'Change in the annual rate of gonorrhea\n(per 100,000 MSM population)'),
        health_measure='u',
        x_range=x_range,
        y_range=y_range,
        cost_multiplier=100000,
        effect_multiplier=SIM_DURATION,
        file_name=fig_file_name,
        fig_size=(5, 5))


