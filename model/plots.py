import apacepy.analysis.scenarios as scen
import apacepy.analysis.trajectories as traj
import apacepy.analysis.visualize_scenarios as vis
from deampy.in_out_functions import write_csv

from definitions import RestProfile, SympStat, REST_PROFILES, ConvertSympAndResitAndAntiBio, \
    SPE_VALUES, SIM_DURATION, ANTIBIOTICS, COVERAGE_VALUES, get_scenario_name
from model import data as D

traj.SUBPLOT_W_SPACE = 0.25
scen.POLY_DEGREES = 1
SCENARIO_COLORS = ['purple', 'blue', 'red', 'green', 'orange', 'brown']


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

    calibration_filename = dir_of_traj_figs+'/(summary) ' + filename

    list_plot_info = [prev, gono_rate, perc_symp]
    list_plot_info.extend(perc_cases_by_rest_profile)
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


def get_rate_percentage_life(scenarios_df, scenario_name):
    """
    :param scenarios_df: (scenario dataframe)
    :param scenario_name: (string) scenario name
    :return: (annual rate of gonorrhea cases, % cases treated with CIP, TET, or CRO, lifespan of CIP, TET, or CRO)
        all calculated after the end of warm-up and formulated as estimated and prediction interval
    """
    rate = scenarios_df.get_mean_interval(
        scenario_name=scenario_name,
        outcome_name='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
        interval_type='p', deci=0, multiplier=100000, form=',')

    proportion = []
    proportion.append(scenarios_df.get_mean_interval(
        scenario_name=scenario_name,
        outcome_name='Time-averaged proportion of cases treated successfully with CIP, TET, or CRO '
                     '(average incidence after epidemic warm-up)',
        interval_type='p', deci=1, form='%'))

    for ab in ANTIBIOTICS:
        proportion.append(scenarios_df.get_mean_interval(
            scenario_name=scenario_name,
            outcome_name='Time-averaged proportion of cases treated successfully with {} '
                         '(average incidence after epidemic warm-up)'.format(ab),
            interval_type='p', deci=1, form='%'))

    life = []
    life.append(scenarios_df.get_mean_interval(
        scenario_name=scenario_name,
        outcome_name='Time-averaged proportion of cases treated successfully with CIP, TET, or CRO '
                     '(average incidence after epidemic warm-up)',
        interval_type='p', deci=1, multiplier=SIM_DURATION))

    for ab in ANTIBIOTICS:
        life.append(
            scenarios_df.get_mean_interval(
                scenario_name=scenario_name,
                outcome_name='Time-averaged proportion of cases treated successfully with {} '
                             '(average incidence after epidemic warm-up)'.format(ab),
                interval_type='p', deci=1, multiplier=SIM_DURATION)
        )

    return rate, proportion, life


def get_perc_change_rate_life(scenarios_df, scenario_name_base, scenario_name_new):

    perc_change_rate = scenarios_df.get_relative_diff_mean_interval(
        scenario_name_base=scenario_name_base,
        scenario_names=scenario_name_new,
        outcome_name='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
        deci=1, form='%')

    perc_change_life = []
    perc_change_life.append(scenarios_df.get_relative_diff_mean_interval(
        scenario_name_base=scenario_name_base,
        scenario_names=scenario_name_new,
        outcome_name='Time-averaged proportion of cases treated successfully with CIP, TET, or CRO '
                     '(average incidence after epidemic warm-up)',
        deci=1, form='%'))

    for ab in ANTIBIOTICS:
        try:
            perc_change_life.append(scenarios_df.get_relative_diff_mean_interval(
                scenario_name_base=scenario_name_base,
                scenario_names=scenario_name_new,
                outcome_name='Time-averaged proportion of cases treated successfully with {} '
                             '(average incidence after epidemic warm-up)'.format(ab),
                deci=1, form='%')
            )
        except:
            perc_change_life.append('')

    return perc_change_rate, perc_change_life


def print_rate_percentage_life(scenarios_df, scenario_name):
    """
       :param scenarios_df: (scenario dataframe)
       :param scenario_name: (string) scenario name
       :return: print
                    annual rate of gonorrhea cases,
                    % cases treated with CIP, TET, or CRO,
                    lifespan of CIP, TET, or CRO)
           all calculated after the end of warm-up and formulated as estimated and confidence interval
       """
    rate, prop, life = get_rate_percentage_life(scenarios_df=scenarios_df, scenario_name=scenario_name)
    print('\nScenario: ', scenario_name)
    print('\tAnnual rate of gonorrhea cases:\t\t', rate)
    print('\tTime-averaged proportion of cases treated with CIP, TET, or CRO:\t', prop)
    print('\tEffective lifespan of CIP, TET, and CFX\t\t', life)


def print_change_rate_percentage_life(scenarios_df, scenario_name_base, scenario_name_new):
    """
       :param scenarios_df: (scenario dataframe)
       :param scenario_name_base: (string) name of the base scenario
       :param scenario_name_new: (string) name of the new scenario
       :return: print
                    change in annual rate of gonorrhea cases,
                    change in % cases treated with CIP, TET, or CRO,
                    change in lifespan of CIP, TET, or CRO)
           all calculated after the end of warm-up and formulated as estimated and confidence interval
       """

    rate = scenarios_df.get_relative_diff_mean_interval(
        scenario_name_base=scenario_name_base,
        scenario_names=scenario_name_new,
        outcome_name='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
        deci=1, form='%')
    life = scenarios_df.get_relative_diff_mean_interval(
        scenario_name_base=scenario_name_base,
        scenario_names=scenario_name_new,
        outcome_name='Time-averaged proportion of cases treated with CIP, TET, or CRO '
                     '(average incidence after epidemic warm-up)',
        deci=1, form='%')

    print('\nChange due to {} compared to {}'.format(scenario_name_new, scenario_name_base))
    print('\tAnnual rate of gonorrhea cases:\t\t', rate)
    print('\tEffective lifespan of CIP, TET, and CFX:\t\t', life)


def export_performance_of_scenarios(if_m_available_for_1st_tx, coverage_values,
                                    simulation_duration=None, calibration_seed=None):
    """
    export performance of different scenarios of test characteristics
    :param if_m_available_for_1st_tx: (bool) if m is available for first-line therapy
    :param simulation_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    :param coverage_values: (float) specific coverage value for the test
    :return: saves a csv file with the following columns for each scenario:
        (rate of gonorrhea cases, % cases treated with CIP, TET, or CRO, lifespan of CIP, TET, or CRO,
        % increase in cases with respect to the status quo,
        % increase in lifespan of CIP, TET, or CRO with respect to the status quo)
    """

    scenario_name = get_scenario_name(
        if_m_available=if_m_available_for_1st_tx,
        sim_duration=simulation_duration,
        calibration_seed=calibration_seed
    )

    csv_file_scenarios = 'outputs/{}/scenarios/simulated_scenarios.csv'.format(scenario_name)
    csv_file_summary = 'outputs/{}/scenarios/performance_summary.csv'.format(scenario_name)

    # read scenarios into a dataframe
    scenarios_df = scen.ScenarioDataFrame(csv_file_name=csv_file_scenarios)

    # rows
    rows = [
        ['Scenario name',
         'Rate of gonorrhea cases',
         '% cases successfully treated with CIP, TET, or CRO',
         '% cases successfully treated with CIP',
         '% cases successfully treated with TET',
         '% cases successfully treated with CRO',
         'Effective lifespan of CIP, TET, and CFX',
         'delta - Rate of gonorrhea cases',
         'delta - Effective lifespan of CIP, TET, and CFX',
         'delta - Effective lifespan of CIP',
         'delta - Effective lifespan of TET',
         'delta - Effective lifespan of CFX',
         ]
    ]

    # status quo
    scenario_name_base = 'Status quo (no rapid test)'
    rate, prob_success, eff_life = get_rate_percentage_life(
        scenarios_df=scenarios_df, scenario_name=scenario_name_base)
    rows.append([scenario_name_base, rate] + prob_success + eff_life + ['', ''])

    #
    for test_coverage in coverage_values:
        # scenario name
        scenario_name = '(p=0.750, q=0.975, c={:.3f})'.format(test_coverage)

        # get rate, percentage treated with 1st-line drugs, and lifespan of 1st-line drugs
        rate, prob_success, eff_life = get_rate_percentage_life(
            scenarios_df=scenarios_df, scenario_name=scenario_name)

        # get % change in rate, % change in lifespan
        perc_change_rate, perc_change_life = get_perc_change_rate_life(
            scenarios_df=scenarios_df, scenario_name_base=scenario_name_base, scenario_name_new=scenario_name)

        # append row
        rows.append([scenario_name, rate] + prob_success + eff_life + [perc_change_rate] + perc_change_life)

    # export to csv
    write_csv(rows=rows, file_name=csv_file_summary)


def get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx, test_coverage=None, sim_duration=None, calibration_seed=None):
    """
    :param if_m_available_for_1st_tx: (bool) if M is available for first-line therapy
    :param test_coverage: (float) coverage of rapid text
    :return: (tuple) name of csv file containing the summary of simulated scenarios,
                     name of figure to save the results as
    :param sim_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    """

    scenario_name = get_scenario_name(
        if_m_available=if_m_available_for_1st_tx,
        sim_duration=sim_duration,
        calibration_seed=calibration_seed)

    if test_coverage is not None:
        fig_file_name = 'figures/SA/{}-coverage {:.2f}.png'.format(scenario_name, test_coverage)
    else:
        fig_file_name = 'figures/SA/{}-coverage.png'.format(scenario_name)

    csv_file_name = 'outputs/{}/scenarios/simulated_scenarios.csv'.format(scenario_name)

    return csv_file_name, fig_file_name


def get_scenarios_with_spec_cov(scenarios_df, i, spec, test_coverage):
    """
    :param scenarios_df: dataframe of simulated scenarios
    :param i: (int) index of the specificity (to get the color code)
    :param spec: (float) specificity
    :param test_coverage: (float) test coverage
    :return: set of scenarios varying sensitivity for the specified specificity and test coverage
    """

    return scen.SetOfScenarios(
        name='Specificity = {:.3f}'.format(spec),
        scenario_df=scenarios_df,
        color=SCENARIO_COLORS[i],
        marker='o',
        conditions_on_variables=[
            scen.ConditionOnVariable(var_name='sensitivity', if_included_in_label=True, label_format='{:.2f}'),
            scen.ConditionOnVariable(var_name='specificity', values=[spec]),
            scen.ConditionOnVariable(var_name='rapid test coverage', values=[test_coverage])
        ],
        if_find_frontier=False,
        if_show_fitted_curve=False,
        labels_shift_x=0.01,
        labels_shift_y=0.01)


def plot_scenarios(csv_file_name, fig_file_name, test_coverage, x_range, y_range, print_all_scenarios=False):
    """ plots the cost-effectiveness figure for a specific text coverage
    :param csv_file_name: (string) csv filename where the summary of simulated scenarios are located
    :param fig_file_name: (string) filename of the figure to save the results as
    :param test_coverage: (float) specific value of test coverage
    :param x_range: range of x-axis
    :param y_range: range of y-axis
    :param print_all_scenarios: (bool) set True to print the performance of different scenarios
    """

    # read scenarios into a dataframe
    scenarios_df = scen.ScenarioDataFrame(csv_file_name=csv_file_name)

    # get a specific outcome from a specific scenario
    if print_all_scenarios:
        print('\nScenario name: ', '\tRate of gonorrhea cases | Effective lifespan of CIP, TET, and CFX')
        for name in scenarios_df.scenarios:
            rate, prop, life = get_rate_percentage_life(scenarios_df=scenarios_df, scenario_name=name)
            print('{}: \t{} | {}'.format(name, rate, life))

    # plot CEA
    scen.ERROR_BAR_ALPHA = 0.2

    # sets of scenarios to display on the cost-effectiveness plain
    list_of_scenario_sets = []

    for i, spec in enumerate(SPE_VALUES):
        list_of_scenario_sets.append(
            get_scenarios_with_spec_cov(scenarios_df=scenarios_df, i=i, spec=spec, test_coverage=test_coverage)
        )

    vis.plot_sets_of_scenarios(
        list_of_scenario_sets=list_of_scenario_sets,
        name_of_base_scenario='Status quo (no rapid test)',
        list_if_remove_base_scenario=[True] * len(SPE_VALUES),
        effect_outcome='Time-averaged proportion of cases treated successfully with CIP, TET, or CRO '
                       '(average incidence after epidemic warm-up)',
        cost_outcome='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
        labels=('Change in the effective lifespan of\nCIP, TET, and CRO (years)',
                'Change in the annual rate of gonorrhea\n(per 100,000 MSM population)'),
        health_measure='u',
        x_range=x_range,
        y_range=y_range,
        cost_multiplier=100000,
        effect_multiplier=SIM_DURATION,
        file_name=fig_file_name,
        fig_size=(5, 5))


def plot_scenario_sa(csv_file_name, fig_file_name, x_range=None, y_range=None, l_b_r_t=None, fig_size=None):

    # read scenarios into a dataframe
    scenarios_df = scen.ScenarioDataFrame(csv_file_name=csv_file_name)

    list_list_series = []
    list_of_titles = []
    for c in COVERAGE_VALUES:
        list_of_titles.append('Test Coverage {:.0%}'.format(c))
        list_of_scenario_sets = []
        for i, spec in enumerate(SPE_VALUES):
            list_of_scenario_sets.append(
                get_scenarios_with_spec_cov(scenarios_df=scenarios_df, i=i, spec=spec, test_coverage=c)
            )
        list_list_series.append(list_of_scenario_sets)

    vis.multi_plot_series(
        list_list_series=list_list_series,
        list_of_titles=list_of_titles,
        name_of_base_scenario='Status quo (no rapid test)',
        list_if_remove_base_scenario=[True] * len(SPE_VALUES),
        effect_outcome='Time-averaged proportion of cases treated successfully with CIP, TET, or CRO '
                       '(average incidence after epidemic warm-up)',
        cost_outcome='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
        x_range=x_range,
        y_range=y_range,
        labels=('Change in the effective lifespan of ciprofloxacin, tetracycline, and ceftriaxone (years)',
                'Change in the annual rate of incident gonorrhea cases\n(per 100,000 MSM population)'),
        cost_multiplier=100000,
        effect_multiplier=SIM_DURATION,
        fig_size=fig_size,
        l_b_r_t=l_b_r_t,
        file_name=fig_file_name
    )
