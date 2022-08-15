import apacepy.analysis.scenarios as scen
import apacepy.analysis.trajectories as traj
import apacepy.analysis.visualize_scenarios as vis
from definitions import RestProfile, SympStat, REST_PROFILES, ConvertSympAndResitAndAntiBio, \
    SIM_DURATION, ANTIBIOTICS, COVERAGE_VALUES, COLOR_VARYING_COVERAGE, \
    CIP_SPEC_VALUES, TET_SPEC_VALUES, COLOR_BY_SPEC, SINGLE_FIG_SIZE, \
    EFFECT_OUTCOME, COST_OUTCOME, EFFECT_COST_LABELS, EFFECT_COST_LABELS_NO_LINE_BREAK
from model import data as D
from model.scenario_and_sensitivity_analyses import get_rate_percentage_life, get_sa_scenarios_with_specific_spec_coverage_ab, \
    get_sa_scenarios_varying_coverage


traj.SUBPLOT_W_SPACE = 0.25
scen.POLY_DEGREES = 1


def plot_trajectories(prev_multiplier=52, incd_multiplier=1,
                      obs_prev_multiplier=1, obs_incd_multiplier=1,
                      dir_of_traj_files='outputs/trajectories',
                      dir_of_traj_figs='figures/trajs',
                      filename='trajectories.png'):
    """
    :param prev_multiplier: (int) to multiply the simulation time to convert it to year, week, or day.
    :param incd_multiplier: (int) to multiply the simulation period to covert it to year, week, or day.
    :param obs_prev_multiplier: (int) to multiply the prevalence survey time to convert it to year, week, or day.
    :param obs_incd_multiplier: (int) to multiply the incidence survey period to covert it to year, week, or day.
    :param dir_of_traj_files: (string) directory where simulated trajectories are located
    :param dir_of_traj_figs: (string) directory where figures should be saved to
    :param filename: (string) filename to save the trajectories as
    """

    sim_outcomes = traj.SimOutcomeTrajectories(csv_directory=dir_of_traj_files)

    # defaults
    traj.TIME_0 = 0  # 2014
    traj.X_RANGE = (0, SIM_DURATION+1)
    traj.X_TICKS = [traj.TIME_0, 5]  # x-axis ticks (min at 0 with interval of 5)
    traj.X_LABEL = 'Year'  # x-axis label
    traj.TRAJ_TRANSPARENCY = 0.75 if 'onetraj' in filename else 0.25

    # plot information
    S = traj.TrajPlotInfo(outcome_name='In: S',
                          title='Susceptible',
                          x_multiplier=prev_multiplier,
                          y_range=(0, 1500000))
    pop = traj.TrajPlotInfo(outcome_name='Population size',
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
            Is.append(traj.TrajPlotInfo(outcome_name='In: I ' + str_symp_susp,
                                        title='I ' + str_symp_susp,
                                        x_multiplier=prev_multiplier))
            # Fs: infectious compartments after treatment failure
            Fs.append(traj.TrajPlotInfo(outcome_name='In: F ' + str_symp_susp,
                                        title='F ' + str_symp_susp,
                                        x_multiplier=prev_multiplier))
            # increment i
            i += 1

    sim_outcomes.plot_multi_panel(n_rows=4, n_cols=4,
                                  list_plot_info=Is,
                                  figure_size=(7, 7),
                                  file_name=dir_of_traj_figs+'/(valid-Is) ' + filename)
    sim_outcomes.plot_multi_panel(n_rows=4, n_cols=4,
                                  list_plot_info=Fs,
                                  figure_size=(7, 7),
                                  file_name=dir_of_traj_figs+'/(valid-Fs) ' + filename)

    # sim_outcomes.plot_multi_panel(n_rows=1, n_cols=2,
    #                               list_plot_info=[S, pop],
    #                               figure_size=(2*2.2, 1*2.2),
    #                               file_name='figures/(valid-S) ' + filename)

    # ------------- Calibration Figure ---------------
    prev = traj.TrajPlotInfo(outcome_name='Prevalence',
                             title='Prevalence (%)',
                             x_multiplier=obs_prev_multiplier, y_multiplier=100,
                             y_range=(0, 10),
                             calibration_info=traj.CalibrationTargetPlotInfo(
                                 rows_of_data=D.Prevalence)
                             )
    gono_rate = traj.TrajPlotInfo(outcome_name='Rate of gonorrhea cases',
                                  title='Rate of gonorrhea cases\n(Per 100,000 MSM population)',
                                  x_multiplier=obs_incd_multiplier, y_multiplier=100000,
                                  y_range=(0, 12500),
                                  calibration_info=traj.CalibrationTargetPlotInfo(
                                      rows_of_data=D.GonorrheaRate)
                                  )

    perc_symp = traj.TrajPlotInfo(outcome_name='Proportion of cases symptomatic',
                                  title='Percent gonorrhea cases\nthat are symptomatic (%)',
                                  x_multiplier=obs_incd_multiplier,
                                  y_multiplier=100, y_range=(0, 100),
                                  calibration_info=traj.CalibrationTargetPlotInfo(
                                      rows_of_data=D.PercSymptomatic)
                                  )

    perc_cases_by_rest_profile = []
    for p in range(len(REST_PROFILES)):
        perc_cases_by_rest_profile.append(traj.TrajPlotInfo(
            outcome_name='Proportion of cases resistant to '+REST_PROFILES[p],
            title='Cases with profile\n{} (%)'.format(REST_PROFILES[p]),
            x_multiplier=obs_incd_multiplier,
            y_multiplier=100, y_range=(0, 100),
            calibration_info=traj.CalibrationTargetPlotInfo(
                rows_of_data=D.PercResistProfile[REST_PROFILES[p]]))
        )

    perc_cases_rest_cro = traj.TrajPlotInfo(
        outcome_name='Proportion of cases CRO-NS',
        title='Cases resistant to CRO',
        x_multiplier=obs_incd_multiplier,
        y_multiplier=100, y_range=(0, 100))

    calibration_filename = dir_of_traj_figs+'/(summary) ' + filename

    list_plot_info = [prev, gono_rate, perc_symp] + perc_cases_by_rest_profile + [perc_cases_rest_cro]

    sim_outcomes.plot_multi_panel(n_rows=4, n_cols=3,
                                  list_plot_info=list_plot_info,
                                  figure_size=(3*2.2, 4*2.3), show_subplot_labels=True,
                                  file_name=calibration_filename)

    # ------------- Successful treatment with different antibiotics ---------------
    Txs = []
    for a in ANTIBIOTICS:
        Txs.append(traj.TrajPlotInfo(outcome_name='To: Counting success tx with ' + a,
                                     title='Successful Tx-' + a,
                                     x_multiplier=incd_multiplier))
    Txs.append(traj.TrajPlotInfo(outcome_name='To: 1st-Tx with M',
                                 title='Successful 1st-Tx-M',
                                 x_multiplier=incd_multiplier))
    Txs.append(traj.TrajPlotInfo(outcome_name='To: Tx with M',
                                 title='Successful Tx-M',
                                 x_multiplier=incd_multiplier))
    sim_outcomes.plot_multi_panel(n_rows=2, n_cols=3,
                                  list_plot_info=Txs,
                                  figure_size=(3*2.2, 2*2.2), show_subplot_labels=True,
                                  file_name=dir_of_traj_figs+'/(valid-Tx) ' + filename)


def plot_sa_for_varying_coverage(csv_file_name, sim_duration, fig_file_name, x_range, y_range):
    """ plots the cost-effectiveness figure for varying test coverage under unknown sensitivity and specificity values
    :param csv_file_name: (string) csv filename where the summary of simulated scenarios are located
    :param sim_duration: (float) simulation duration
    :param fig_file_name: (string) filename of the figure to save the results as
    :param x_range: range of x-axis
    :param y_range: range of y-axis
    """

    # read scenarios into a dataframe
    scenarios_df = scen.ScenarioDataFrame(csv_file_name=csv_file_name)

    # plot CEA
    scen.ERROR_BAR_ALPHA = 0.2

    # sets of scenarios to display on the cost-effectiveness plane
    list_of_scenario_sets = [get_sa_scenarios_varying_coverage(
        scenarios_df=scenarios_df, color=COLOR_VARYING_COVERAGE)]

    vis.plot_sets_of_scenarios(
        list_of_scenario_sets=list_of_scenario_sets,
        name_of_base_scenario='Status quo (no rapid test)',
        list_if_remove_base_scenario=[True],
        effect_outcome=EFFECT_OUTCOME,
        cost_outcome=COST_OUTCOME,
        labels=EFFECT_COST_LABELS,
        health_measure='u',
        x_range=x_range,
        y_range=y_range,
        cost_multiplier=100000,
        effect_multiplier=sim_duration,
        file_name=fig_file_name,
        fig_size=SINGLE_FIG_SIZE)


def plot_sa_for_specific_ab_and_coverage(csv_file_name, fig_file_name, ab, test_coverage,
                                         x_range, y_range, print_all_scenarios=False):
    """ plots the cost-effectiveness figure for a specific antibiotic and test coverage
    :param csv_file_name: (string) csv filename where the summary of simulated scenarios are located
    :param fig_file_name: (string) filename of the figure to save the results as
    :param ab: (string) 'CIP' or 'TET'
    :param test_coverage: (float) specific value of test coverage
    :param x_range: range of x-axis
    :param y_range: range of y-axis
    :param print_all_scenarios: (bool) set True to print the performance of different scenarios
    """

    # read scenarios into a dataframe
    scenarios_df = scen.ScenarioDataFrame(csv_file_name=csv_file_name)

    # get a specific outcome from a specific scenario
    if print_all_scenarios:
        print('\nScenario name: ', '\tRate of gonorrhea cases | Effective lifespan of CIP, TET, and CRO')
        for name in scenarios_df.scenarios:
            rate, prop, life = get_rate_percentage_life(scenarios_df=scenarios_df, scenario_name=name)
            print('{}: \t{} | {}'.format(name, rate, life))

    # plot CEA
    scen.ERROR_BAR_ALPHA = 0.2

    # sets of scenarios to display on the cost-effectiveness plane
    list_of_scenario_sets = []

    if ab == 'CIP':
        spec_values = CIP_SPEC_VALUES
    elif ab == 'TET':
        spec_values = TET_SPEC_VALUES
    else:
        raise ValueError('Invalid value for antibiotic.')

    for i, spec in enumerate(spec_values):
        list_of_scenario_sets.append(
            get_sa_scenarios_with_specific_spec_coverage_ab(
                scenarios_df=scenarios_df, spec=spec, test_coverage=test_coverage, ab=ab, color=COLOR_BY_SPEC[i])
        )

    vis.plot_sets_of_scenarios(
        list_of_scenario_sets=list_of_scenario_sets,
        name_of_base_scenario='Status quo (no rapid test)',
        list_if_remove_base_scenario=[True] * len(spec_values),
        effect_outcome=EFFECT_OUTCOME,
        cost_outcome=COST_OUTCOME,
        labels=EFFECT_COST_LABELS,
        health_measure='u',
        x_range=x_range,
        y_range=y_range,
        cost_multiplier=100000,
        effect_multiplier=SIM_DURATION,
        file_name=fig_file_name,
        fig_size=SINGLE_FIG_SIZE)


def plot_sa_for_specific_ab(ab, csv_file_name, fig_file_name,
                            x_range=None, y_range=None, l_b_r_t=None, fig_size=None):

    # read scenarios into a dataframe
    scenarios_df = scen.ScenarioDataFrame(csv_file_name=csv_file_name)

    if ab == 'CIP':
        spec_values = CIP_SPEC_VALUES
    elif ab == 'TET':
        spec_values = TET_SPEC_VALUES
    else:
        raise ValueError('Invalid value for antibiotic.')

    list_list_series = []
    list_of_titles = []
    for c in COVERAGE_VALUES:
        list_of_titles.append('Test Coverage {:.0%}'.format(c))
        list_of_scenario_sets = []
        for i, spec in enumerate(spec_values):
            list_of_scenario_sets.append(
                get_sa_scenarios_with_specific_spec_coverage_ab(
                    scenarios_df=scenarios_df, spec=spec, test_coverage=c, ab=ab, color=COLOR_BY_SPEC[i])
            )
        list_list_series.append(list_of_scenario_sets)

    vis.multi_plot_series(
        list_list_series=list_list_series,
        list_of_titles=list_of_titles,
        name_of_base_scenario='Status quo (no rapid test)',
        list_if_remove_base_scenario=[True] * len(spec_values),
        effect_outcome=EFFECT_OUTCOME,
        cost_outcome=COST_OUTCOME,
        x_range=x_range,
        y_range=y_range,
        labels=EFFECT_COST_LABELS_NO_LINE_BREAK,
        cost_multiplier=100000,
        effect_multiplier=SIM_DURATION,
        fig_size=fig_size,
        l_b_r_t=l_b_r_t,
        file_name=fig_file_name
    )
