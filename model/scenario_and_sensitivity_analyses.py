import apacepy.analysis.scenarios as scen
from deampy.in_out_functions import write_csv
from deampy.sensitivity_analysis import SensitivityAnalysis

from definitions import ROOT_DIR, get_scenario_name, ANTIBIOTICS, get_name_of_scenario_analysis, \
    TRANSMISSION_FACTOR_VALUES, FIG_EXT

SCENARIO_COLORS = ['purple', 'blue', 'red', 'green', 'orange', 'brown']


def get_rate_percentage_life(scenarios_df, scenario_name, sim_duration):
    """
    :param scenarios_df: (scenario dataframe)
    :param scenario_name: (string) scenario name
    :param sim_duration: (float) simulation duration
    :return: (annual rate of gonorrhea cases, % cases treated with CIP, TET, or CRO, lifespan of CIP, TET, or CRO)
        all calculated after the end of warm-up and formulated as estimated and prediction interval
    """
    rate = scenarios_df.get_mean_interval(
        scenario_name=scenario_name,
        outcome_name='Rate of gonorrhea cases (average incidence after epidemic warm-up)',
        interval_type='p', deci=0, multiplier=100000, form=',')

    proportion = [scenarios_df.get_mean_interval(
        scenario_name=scenario_name,
        outcome_name='Time-averaged proportion of cases treated successfully with CIP, TET, or CRO '
                     '(average incidence after epidemic warm-up)',
        interval_type='p', deci=1, form='%')]

    for ab in ANTIBIOTICS:
        proportion.append(scenarios_df.get_mean_interval(
            scenario_name=scenario_name,
            outcome_name='Time-averaged proportion of cases {}-S '
                         '(average incidence after epidemic warm-up)'.format(ab),
            interval_type='p', deci=1, form='%'))

    life = [scenarios_df.get_mean_interval(
        scenario_name=scenario_name,
        outcome_name='Time-averaged proportion of cases treated successfully with CIP, TET, or CRO '
                     '(average incidence after epidemic warm-up)',
        interval_type='p', deci=1, multiplier=sim_duration)]

    for ab in ANTIBIOTICS:
        life.append(
            scenarios_df.get_mean_interval(
                scenario_name=scenario_name,
                outcome_name='Time-averaged proportion of cases {}-S '
                             '(average incidence after epidemic warm-up)'.format(ab),
                interval_type='p', deci=1, multiplier=sim_duration)
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
                outcome_name='Time-averaged proportion of cases {}-S '
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
                                    simulation_duration, calibration_seed=None,
                                    varying_trans_factor=False, if_wider_priors=False):
    """
    export performance of different scenarios of test characteristics
    :param if_m_available_for_1st_tx: (bool) if m is available for first-line therapy
    :param coverage_values: (float) specific coverage value for the test
    :param simulation_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    :param varying_trans_factor: (float) if transmission factor is being varied in this analysis
    :param if_wider_priors: (bool) if model with wider priors should be used
    :return: saves a csv file with the following columns for each scenario:
        (rate of gonorrhea cases, % cases treated with CIP, TET, or CRO, lifespan of CIP, TET, or CRO,
        % increase in cases with respect to the status quo,
        % increase in lifespan of CIP, TET, or CRO with respect to the status quo)
    """

    scenario_name = get_scenario_name(
        if_m_available=if_m_available_for_1st_tx,
        sim_duration=simulation_duration,
        calibration_seed=calibration_seed,
        if_varying_transmission_factor=varying_trans_factor,
        if_wider_priors=if_wider_priors
    )

    csv_file_scenarios = 'outputs/scen-{}/scenarios/simulated_scenarios.csv'.format(scenario_name)
    csv_file_summary = 'outputs/scen-{}/scenarios/performance_summary.csv'.format(scenario_name)

    # read scenarios into a dataframe
    scenarios_df = scen.ScenarioDataFrame(csv_file_name=csv_file_scenarios)

    # rows
    rows = [
        ['Scenario name',
         'Rate of gonorrhea cases',
         '% cases successfully treated with CIP, TET, or CRO',
         '% cases susceptible to CIP',
         '% cases susceptible to TET',
         '% cases susceptible to CRO',
         'Effective lifespan of CIP, TET, and CRO',
         'Effective lifespan of CIP',
         'Effective lifespan of TET',
         'Effective lifespan of CRO',
         'delta - Rate of gonorrhea cases',
         'delta - Effective lifespan of CIP, TET, and CRO',
         'delta - Effective lifespan of CIP',
         'delta - Effective lifespan of TET',
         'delta - Effective lifespan of CRO',
         ]
    ]

    # status quo
    scenario_name_base = 'Status quo (no rapid test)'
    rate, prob_success, eff_life = get_rate_percentage_life(
        scenarios_df=scenarios_df, scenario_name=scenario_name_base, sim_duration=simulation_duration)
    rows.append([scenario_name_base, rate] + prob_success + eff_life + ['', '', '', '', ''])

    # find the range of transmission factors
    trans_factors = [1.0] if not varying_trans_factor else TRANSMISSION_FACTOR_VALUES

    for f in trans_factors:
        for test_coverage in coverage_values:
            # scenario name
            scenario_name = get_name_of_scenario_analysis(
                cip_sens=None, cip_spec=None, tet_sens=None, tet_spec=None,
                coverage=test_coverage, transmission_factor=f)

            # get rate, percentage treated with 1st-line drugs, and lifespan of 1st-line drugs
            rate, prob_success, eff_life = get_rate_percentage_life(
                scenarios_df=scenarios_df, scenario_name=scenario_name, sim_duration=simulation_duration)

            # get % change in rate, % change in lifespan
            perc_change_rate, perc_change_life = get_perc_change_rate_life(
                scenarios_df=scenarios_df, scenario_name_base=scenario_name_base, scenario_name_new=scenario_name)

            # append row
            rows.append([scenario_name, rate] + prob_success + eff_life + [perc_change_rate] + perc_change_life)

    # export to csv
    write_csv(rows=rows, file_name=csv_file_summary)


def get_scenarios_csv_filename_and_fig_filename(
        if_m_available_for_1st_tx, ab=None, test_coverage=None,
        sim_duration=None, calibration_seed=None, varying_trans_factor=False, if_wider_priors=False):
    """
    :param if_m_available_for_1st_tx: (bool) if M is available for first-line therapy
    :param ab: (string) 'CIP' or 'TET'
    :param test_coverage: (float) coverage of rapid test
    :param sim_duration: (float) simulation duration (for sensitivity analysis)
    :param calibration_seed: (int) calibration seed (for sensitivity analysis)
    :param varying_trans_factor: (float) if varying the transmission factor
    :param if_wider_priors: (bool) if using the model with wider prior distributions
    :return: (tuple) name of csv file containing the summary of simulated scenarios,
                 name of figure to save the results as
    """

    scenario_name = get_scenario_name(
        if_m_available=if_m_available_for_1st_tx,
        sim_duration=sim_duration,
        calibration_seed=calibration_seed,
        if_varying_transmission_factor=varying_trans_factor,
        if_wider_priors=if_wider_priors)

    if test_coverage is not None:
        fig_file_name = 'figures/SA/{}-{}-coverage {:.2f}.{}'.format(scenario_name, ab, test_coverage, FIG_EXT)
    else:
        if ab is None:
            fig_file_name = 'figures/SA/{}.{}'.format(scenario_name, FIG_EXT)
        else:
            fig_file_name = 'figures/SA/{}-{}-coverage.{}'.format(scenario_name, ab, FIG_EXT)

    csv_file_name = 'outputs/scen-{}/scenarios/simulated_scenarios.csv'.format(scenario_name)

    return csv_file_name, fig_file_name


def get_sa_scenarios_varying_coverage(scenarios_df, color, transmission_factor=1.0):
    """
    :param scenarios_df: dataframe of simulated scenarios
    :param color: (string) color of set
    :param transmission_factor: (float) transmission factor
    :return: set of scenarios varying sensitivity for the specified specificity and test coverage
    """

    return scen.SetOfScenarios(
        name='Transmission factor {:.2f}'.format(transmission_factor),
        scenario_df=scenarios_df,
        color=color,
        marker='o',
        conditions_on_variables=[
            scen.ConditionOnVariable(var_name='CIP-sens', values=[None, '']),
            scen.ConditionOnVariable(var_name='CIP-spec', values=[None, '']),
            scen.ConditionOnVariable(var_name='TET-sens', values=[None, '']),
            scen.ConditionOnVariable(var_name='TET-spec', values=[None, '']),
            scen.ConditionOnVariable(var_name='rapid test coverage', if_included_in_label=True, label_format='{:.0%}'),
            scen.ConditionOnVariable(var_name='transmission factor', values=transmission_factor)
        ],
        if_find_frontier=False,
        if_show_fitted_curve=False,
        labels_shift_x=0.035,
        labels_shift_y=-0.01)


def get_sa_scenarios_with_specific_spec_coverage_ab(
        scenarios_df, spec, test_coverage, ab, color, include_sens_labels):
    """
    :param scenarios_df: dataframe of simulated scenarios
    :param spec: (float) specificity
    :param test_coverage: (float) test coverage
    :param ab: (string): 'CIP' or 'TET'
    :param color: (string) color of set
    :param include_sens_labels: (bool) set to True to show the sensitivity values on the figure
    :return: set of scenarios varying sensitivity for the specified specificity and test coverage
    """

    if ab == 'CIP':
        other_ab = 'TET'
    elif ab == 'TET':
        other_ab = 'CIP'
    else:
        raise ValueError('Invalid value for antibiotic.')

    return scen.SetOfScenarios(
        name='Specificity = {:.2f}'.format(spec),
        scenario_df=scenarios_df,
        color=color,
        marker='o',
        conditions_on_variables=[
            scen.ConditionOnVariable(var_name='{}-sens'.format(other_ab), values=[None, '']),
            scen.ConditionOnVariable(var_name='{}-spec'.format(other_ab), values=[None, '']),
            scen.ConditionOnVariable(var_name='{}-sens'.format(ab),
                                     if_included_in_label=include_sens_labels, label_format='{:.2f}'),
            scen.ConditionOnVariable(var_name='{}-spec'.format(ab), values=[spec]),
            scen.ConditionOnVariable(var_name='rapid test coverage', values=[test_coverage])
        ],
        if_find_frontier=False,
        if_show_fitted_curve=False,
        labels_shift_x=0.015,
        labels_shift_y=0.01)


def print_corr(scenarios, dict_select_params, outcome, outcome_label):

    # find the increase in effective lifespan
    increase = scenarios.scenarios['(p=None, None), q=(None, None), c=0.750)'].outcomes[outcome] - \
               scenarios.scenarios['Status quo (no rapid test)'].outcomes[outcome]
    print('Increase in {}:'.format(outcome_label), sum(increase) / len(increase))

    # calculate the partial correlation coefficients of all parameters
    sa = SensitivityAnalysis(dic_parameter_values=dict_select_params,
                             output_values=increase)

    sa.export_to_csv(corrs=['r', 'p', 'pr'], decimal=4,
                     file_name=ROOT_DIR + '/analysis/outputs/partial_corr_{}.csv'.format(outcome_label))
