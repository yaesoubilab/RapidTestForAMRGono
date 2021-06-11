from SimPy.Parameters import Constant
from apace.CalibrationSupport import FeasibleConditions
from apace.Compartment import Compartment, ChanceNode
from apace.Control import InterventionAffectingEvents, ConditionBasedDecisionRule
from apace.Event import EpiIndepEvent, EpiDepEvent
from apace.FeaturesAndConditions import FeatureSurveillance, FeatureIntervention, \
    ConditionOnFeatures, ConditionOnConditions, ConditionAlwaysFalse

from apace.TimeSeries import SumPrevalence, SumIncidence, RatioTimeSeries
from model.ModelParameters import Parameters
from definitions import SuspProfile, AB, SympStat, SUSP_PROFILES, ComboSympAndSusp


def build_model(model):

    # model settings
    sets = model.settings
    indexer = ComboSympAndSusp(n_symp_stats=len(SympStat), n_susp_profiles=len(SuspProfile))

    # model parameters
    params = Parameters(model_sets=sets)

    # ------------- model compartments ---------------
    Is = [None] * indexer.length
    Ts = [None] *
    RTs = [None] * indexer.length

    # model compartments
    S = Compartment(name='S', size_par=params.sizeS,
                    susceptibility_params=[Constant(value=1)])

    for s in range(len(SympStat)):
        for p in range(len(SuspProfile)):
            i = indexer.get_row_index(symp_state=s, susp_profile=p)
            infectivity_params = [Constant(value=0)] * indexer.nSuspProfiles
            infectivity_params[p] = params.infectivityBySuspProfile[p]
            Is[i] = Compartment(name='I'+indexer.get_str_susp_profile(symp_state=s, susp_profile=p),
                                size_par=params.sizeIBySympAndSusp[i],
                                if_empty_to_eradicate=True,
                                infectivity_params=infectivity_params)

    I_Sus_Sym = I(sus_profile='Sus', symptom='Sym', params=params)
    I_Sus_Asym = I(sus_profile='Sus', symptom='Asym', params=params)
    I_Res_Sym = I(sus_profile='Res', symptom='Sym', params=params)
    I_Res_Asym = I(sus_profile='Res', symptom='Asym', params=params)

    W_I_Sus_Sym = Compartment(name='Waiting 1st | '+I_Sus_Sym.name, num_of_pathogens=2)
    W_I_Sus_Asym = Compartment(name='Waiting 1st | '+I_Sus_Asym.name, num_of_pathogens=2)
    W_I_Res_Sym = Compartment(name='Waiting 1st | '+I_Res_Sym.name, num_of_pathogens=2)
    W_I_Res_Asym = Compartment(name='Waiting 1st | '+I_Res_Asym.name, num_of_pathogens=2)
    W2_I_Res_Sym = Compartment(name='Waiting 2nd | ' + I_Res_Sym.name, num_of_pathogens=2)

    if_Sym_Sus = ChanceNode(name='If symptomatic | Sus',
                            destination_compartments=[I_Sus_Sym, I_Sus_Asym],
                            probability_params=params.dictOfParams['Prob symptomatic'])
    if_sym_inf_res = ChanceNode(name='If symptomatic | Res',
                                destination_compartments=[I_Res_Sym, I_Res_Asym],
                                probability_params=params.dictOfParams['Prob symptomatic'])
    if_sym_emerge_res = ChanceNode(name='If symptomatic after emergence of resistance',
                                   destination_compartments=[W2_I_Res_Sym, I_Res_Asym],
                                   probability_params=params.dictOfParams['Prob symptomatic'])

    if_res_after_tx = ChanceNode(name='If resistance after Tx-A',
                                 destination_compartments=[if_sym_emerge_res, S],
                                 probability_params=params.dictOfParams['Prob of resistance'])

    # ------------- compartment histories ---------------
    # set up prevalence, incidence, and cumulative incidence to collect
    S.setup_history(collect_prev=True)
    I_Sus_Sym.setup_history(collect_prev=True, collect_incd=True)
    I_Sus_Asym.setup_history(collect_prev=True, collect_incd=True)
    I_Res_Sym.setup_history(collect_prev=True, collect_incd=True)
    I_Res_Asym.setup_history(collect_prev=True, collect_incd=True)
    W2_I_Res_Sym.setup_history(collect_prev=True)

    # ------------- summation statistics ---------------
    # population size
    all_comparts = [S,
                    I_Sus_Sym, I_Sus_Asym, I_Res_Sym, I_Res_Asym,
                    W_I_Sus_Sym, W_I_Sus_Asym, W_I_Res_Sym, W_I_Res_Asym,
                    W2_I_Res_Sym]
    pop_size = SumPrevalence(name='Population size', compartments=all_comparts)

    # number infected
    n_infected = SumPrevalence(name='Number infected',
                               compartments=[I_Sus_Sym, I_Sus_Asym, I_Res_Sym, I_Res_Asym,
                                             W_I_Sus_Sym, W_I_Sus_Asym, W_I_Res_Sym, W_I_Res_Asym,
                                             W2_I_Res_Sym])
    prevalence = RatioTimeSeries(name='Prevalence',
                                 numerator_sum_time_series=n_infected,
                                 denominator_sum_time_series=pop_size,
                                 if_surveyed=True)

    # rate of gonorrhea cases
    new_cases = SumIncidence(name='Number of gonorrhea cases',
                             compartments=[W_I_Sus_Sym, W_I_Sus_Asym, W_I_Res_Sym, W_I_Res_Asym])
    gono_rate = RatioTimeSeries(name='Rate of gonorrhea cases',
                                numerator_sum_time_series=new_cases,
                                denominator_sum_time_series=pop_size,
                                if_surveyed=True)

    # % cases symptomatic
    cases_symptomatic = SumIncidence(name='Number of cases symptomatic',
                                     compartments=[W_I_Sus_Sym, W_I_Res_Sym])
    percent_cases_symptomatic = RatioTimeSeries(name='Proportion of cases symptomatic',
                                                numerator_sum_time_series=cases_symptomatic,
                                                denominator_sum_time_series=new_cases,
                                                if_surveyed=True)

    # cases resistant to drug A
    cases_rest_A = SumIncidence(name='Number of cases resistant',
                                compartments=[W_I_Res_Sym, W_I_Res_Asym])
    percent_cases_rest_A = RatioTimeSeries(
        name='Proportion of cases resistant',
        numerator_sum_time_series=cases_rest_A,
        denominator_sum_time_series=new_cases,
        if_surveyed=True,
        survey_size_param=params.dictOfParams['Survey size (over observation periods)'])

    # ------------- calibration targets ---------------
    if sets.calcLikelihood:
        # prevalence
        prevalence.add_feasible_conditions(feasible_conditions=FeasibleConditions(
            feasible_min=0, feasible_max=0.1))
        prevalence.add_calibration_targets(ratios=sets.prevMean,
                                           survey_sizes=sets.prevN)
        # gonorrhea rate
        gono_rate.add_feasible_conditions(feasible_conditions=FeasibleConditions(
            feasible_min=0.04, feasible_max=0.1))
        gono_rate.add_calibration_targets(ratios=sets.gonoRateMean,
                                          survey_sizes=sets.gonoRateN)
        # % cases symptomatic
        percent_cases_symptomatic.add_feasible_conditions(feasible_conditions=FeasibleConditions(
            feasible_min=0.5, feasible_max=1))
        percent_cases_symptomatic.add_calibration_targets(ratios=sets.percSympMean,
                                                          survey_sizes=sets.percSympN)

    # ------------- interventions ---------------
    # interventions (first line therapy with A or M)
    first_line_tx_with_A = InterventionAffectingEvents(name='1st line therapy with Drug A')
    first_line_tx_with_M = InterventionAffectingEvents(name='1st line therapy with Drug M')

    # ------------- features ---------------
    # features
    f_perc_resist_A = FeatureSurveillance(name='Surveyed % of cases resistant',
                                          ratio_time_series_with_surveillance=percent_cases_rest_A)
    f_if_M_ever_switched_on = FeatureIntervention(name='If Drug M ever switched on',
                                                  intervention=first_line_tx_with_M,
                                                  feature_type='if ever switched on')

    # ------------- conditions ---------------
    # conditions
    threshold = 0.1
    A_in_condition = ConditionOnFeatures(name='If % resistant is below threshold',
                                         features=[f_perc_resist_A],
                                         signs=['l'],
                                         thresholds=[threshold])
    A_out_condition = ConditionOnFeatures(name='If % resistant passes threshold',
                                          features=[f_perc_resist_A],
                                          signs=['ge'],
                                          thresholds=[threshold])
    M_is_never_used = ConditionOnFeatures(name='If M is ever used as 1st-line therapy',
                                          features=[f_if_M_ever_switched_on],
                                          signs=['e'],
                                          thresholds=[0])

    turn_on_tx_with_A = ConditionOnConditions(name='', conditions=[A_in_condition, M_is_never_used])

    # ------------- decision rules ---------------
    # add decision rules to interventions
    first_line_tx_with_A.add_decision_rule(
        decision_rule=ConditionBasedDecisionRule(default_switch_value=1,
                                                 condition_to_turn_on=turn_on_tx_with_A,
                                                 condition_to_turn_off=A_out_condition))
    first_line_tx_with_M.add_decision_rule(
        decision_rule=ConditionBasedDecisionRule(default_switch_value=0,
                                                 condition_to_turn_on=A_out_condition,
                                                 condition_to_turn_off=ConditionAlwaysFalse()))

    # ------------- attach epidemic events ---------------
    # attached epidemic events to compartments
    # add events to S
    S.add_events(events=[
        infection(sus_profile='Sus', destination=if_Sym_Sus),
        infection(sus_profile='Res', destination=if_sym_inf_res)
    ])
    # add events to infection compartments
    add_events_to_I(compart=I_Sus_Sym, symptom='Sym', susceptible_compart=S, waiting_compart=W_I_Sus_Sym, params=params)
    add_events_to_I(compart=I_Sus_Asym, symptom='Asym', susceptible_compart=S, waiting_compart=W_I_Sus_Asym, params=params)
    add_events_to_I(compart=I_Res_Sym, symptom='Sym', susceptible_compart=S, waiting_compart=W_I_Res_Sym, params=params)
    add_events_to_I(compart=I_Res_Asym, symptom='Asym', susceptible_compart=S, waiting_compart=W_I_Res_Asym, params=params)
    # add events to "waiting for treatment" compartments
    W_I_Sus_Sym.add_event(
        event=Tx1st(compart=W_I_Sus_Sym, drug='A', destination=if_res_after_tx, intervention=first_line_tx_with_A))
    W_I_Sus_Sym.add_event(
        event=Tx1st(compart=W_I_Sus_Sym, drug='M', destination=S, intervention=first_line_tx_with_M))
    W_I_Sus_Asym.add_event(
        event=Tx1st(compart=W_I_Sus_Asym, drug='A', destination=if_res_after_tx, intervention=first_line_tx_with_A))
    W_I_Sus_Asym.add_event(
        event=Tx1st(compart=W_I_Sus_Asym, drug='M', destination=S, intervention=first_line_tx_with_M))
    W_I_Res_Sym.add_event(
        event=Tx1st(compart=W_I_Res_Sym, drug='A', destination=W2_I_Res_Sym, intervention=first_line_tx_with_A))
    W_I_Res_Sym.add_event(
        event=Tx1st(compart=W_I_Res_Sym, drug='M', destination=S, intervention=first_line_tx_with_M))
    W_I_Res_Asym.add_event(
        event=Tx1st(compart=W_I_Res_Asym, drug='A', destination=I_Res_Asym, intervention=first_line_tx_with_A))
    W_I_Res_Asym.add_event(
        event=Tx1st(compart=W_I_Res_Asym, drug='M', destination=S, intervention=first_line_tx_with_M))

    W2_I_Res_Sym.add_event(
        event=Tx2nd(compart=W2_I_Res_Sym, destination=S, params=params))

    # ------------- population the model ---------------
    # populate the model
    model.populate(compartments=all_comparts,
                   chance_nodes=[if_Sym_Sus, if_sym_inf_res, if_sym_emerge_res, if_res_after_tx],
                   list_of_sum_time_series=[pop_size, n_infected, new_cases, cases_symptomatic, cases_rest_A],
                   list_of_ratio_time_series=[prevalence, gono_rate, percent_cases_symptomatic, percent_cases_rest_A],
                   interventions=[first_line_tx_with_A, first_line_tx_with_M],
                   features=[f_perc_resist_A, f_if_M_ever_switched_on],
                   conditions=[A_in_condition, A_out_condition, M_is_never_used, ],
                   parameters=params)


def I(sus_profile, symptom, params):
    """
    :param sus_profile: (string) 'Sus' for susceptible and 'Res' for resistant
    :param symptom: (string) 'Sym' or 'Asym'
    :param params: (Parameters) model parameters
    :return: an infectious compartment
    """

    if sus_profile == 'Sus':
        infectivity_params = [params.dictOfParams['Infectivity of susceptible strain'],
                              Constant(value=0)]
    elif sus_profile == 'Res':
        infectivity_params = [Constant(value=0),
                              params.dictOfParams['Infectivity of resistant strain']]
    else:
        raise ValueError('Invalid profile')

    return Compartment(name='I | {} | {}'.format(sus_profile, symptom),
                       size_par=params.dictOfParams['Size of I-{}|{}'.format(sus_profile, symptom)],
                       infectivity_params=infectivity_params,
                       if_empty_to_eradicate=True)


def infection(sus_profile, destination):
    """
    :param sus_profile: (string) 'S' for susceptible and 'R' for resistant
    :param destination: (Compartment) destination
    :return: an epidemic dependent event
    """
    pathogen = 0 if sus_profile == 'Sus' else 1

    return EpiDepEvent(
        name='Infection | ' + sus_profile, destination=destination, generating_pathogen=pathogen)


def natural_recovery(compart, destination, params):
    return EpiIndepEvent(
        name='Natural recovery | ' + compart.name,
        rate_param=params.dictOfParams['Rate of natural recovery'],
        destination=destination)


def seeking_treatment(compart, destination, params):
    return EpiIndepEvent(
        name='Seeking treatment | ' + compart.name,
        rate_param=params.dictOfParams['Rate of seeking treatment'],
        destination=destination)


def screened(compart, destination, params):
    return EpiIndepEvent(
        name='Screening | ' + compart.name,
        rate_param=params.dictOfParams['Rate of screening'],
        destination=destination)


def Tx1st(compart, drug, destination, intervention):
    return EpiIndepEvent(
        name='1st Tx with {} | '.format(drug) + compart.name,
        rate_param=Constant(value=1000),
        destination=destination,
        interv_to_activate=intervention)


def Tx2nd(compart, destination, params):
    return EpiIndepEvent(
        name='2nd Tx with M | ' + compart.name,
        rate_param=params.dictOfParams['Rate of seeking retreatment'],
        destination=destination)


def add_events_to_I(compart, symptom, susceptible_compart, waiting_compart, params):

    events = [natural_recovery(compart=compart, params=params, destination=susceptible_compart),
              screened(compart=compart, params=params, destination=waiting_compart)]
    if symptom == 'Sym':
        events.append(seeking_treatment(compart=compart, params=params, destination=waiting_compart))

    compart.add_events(events=events)