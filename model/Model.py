from SimPy.Parameters import Constant
from apace.CalibrationSupport import FeasibleConditions
from apace.Compartment import Compartment, ChanceNode
from apace.Control import InterventionAffectingEvents, ConditionBasedDecisionRule
from apace.Event import EpiIndepEvent, EpiDepEvent
from apace.FeaturesAndConditions import FeatureSurveillance, FeatureIntervention, \
    ConditionOnFeatures, ConditionOnConditions, ConditionAlwaysFalse
from apace.TimeSeries import SumPrevalence, SumIncidence, RatioTimeSeries
from definitions import RestProfile, AB, SympStat, REST_PROFILES, ConvertSympAndSuspAndAntiBio
from model.ModelParameters import Parameters


def build_model(model):

    # model settings
    sets = model.settings
    covert_symp_susp = ConvertSympAndSuspAndAntiBio(
        n_symp_stats=len(SympStat), n_susp_profiles=len(RestProfile))
    covert_symp_susp_antibio = ConvertSympAndSuspAndAntiBio(
        n_symp_stats=len(SympStat), n_susp_profiles=len(RestProfile), n_antibiotics=len(AB))

    # model parameters
    params = Parameters(model_sets=sets)

    # ------------- model compartments ---------------
    Is = [None] * covert_symp_susp.length
    Fs = [None] * covert_symp_susp.length
    ifs_symp_from_S = [None] * len(RestProfile)
    ifs_re_tx = [None] * covert_symp_susp.length
    ifs_symp_from_emerg_res = [None] * len(RestProfile)
    ifs_tx_outcomes = [None] * covert_symp_susp_antibio.length
    ifs_rapid_tests = [None] * covert_symp_susp.length
    Is_symp = []
    if_rest_to = [[]] * len(RestProfile)

    # S
    S = Compartment(name='S', size_par=params.sizeS,
                    susceptibility_params=[Constant(value=1)])

    # Is (infectious) and Fs (infectious after treatment failure)
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            i = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
            # Is
            name = 'I '+covert_symp_susp.get_str_symp_susp(symp_state=s, susp_profile=p)
            infectivity_params = [Constant(value=0)] * covert_symp_susp.nSuspProfiles
            infectivity_params[p] = params.infectivityBySuspProfile[p]
            Is[i] = Compartment(name=name,
                                size_par=params.sizeIBySympAndSusp[i],
                                if_empty_to_eradicate=True,
                                infectivity_params=infectivity_params)
            # if symptomatic
            if s == SympStat.SYMP:
                Is_symp.append(Is[i])

            # Fs: infectious compartments after treatment failure
            name = 'F ' + covert_symp_susp.get_str_symp_susp(symp_state=s, susp_profile=p)
            infectivity_params = [Constant(value=0)] * covert_symp_susp.nSuspProfiles
            infectivity_params[p] = params.infectivityBySuspProfile[p]
            Fs[i] = Compartment(name=name,
                                if_empty_to_eradicate=True,
                                infectivity_params=infectivity_params)

    # if symptomatic after infection
    for p in range(len(RestProfile)):
        name = 'If symptomatic to ' + covert_symp_susp.get_str_susp(susp_profile=p)
        dest_symp = Is[covert_symp_susp.get_row_index(symp_state=SympStat.SYMP.value, susp_profile=p)]
        dest_asym = Is[covert_symp_susp.get_row_index(symp_state=SympStat.ASYM.value, susp_profile=p)]
        ifs_symp_from_S[p] = ChanceNode(name=name,
                                        destination_compartments=[dest_symp, dest_asym],
                                        probability_params=params.probSym)

    # if retreatment after treatment failure
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            name = 'Re-Tx after failure in ' + covert_symp_susp.get_str_symp_susp(symp_state=s, susp_profile=p)
            j = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
            if s == SympStat.SYMP.value:
                dest_yes = Fs[j]
                dest_no = Is[j]
                prob_yes = Constant(1)
            elif s == SympStat.ASYM.value:
                dest_yes = Fs[j]
                dest_no = Is[j]
                prob_yes = Constant(0)
            else:
                raise ValueError('Invalid symptom status.')
            i = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
            ifs_re_tx[i] = ChanceNode(name=name,
                                      destination_compartments=[dest_yes, dest_no],
                                      probability_params=prob_yes)

    # chance nodes for if symptomatic after emergence of resistance during treatment
    for p in range(len(RestProfile)):
        name = 'If symptomatic after moving to {} due to the emergence of resistance'.format(
            covert_symp_susp.get_str_susp(susp_profile=p))
        dest_symp = ifs_re_tx[covert_symp_susp.get_row_index(symp_state=SympStat.SYMP.value, susp_profile=p)]
        dest_asym = ifs_re_tx[covert_symp_susp.get_row_index(symp_state=SympStat.ASYM.value, susp_profile=p)]
        ifs_symp_from_emerg_res[p] = ChanceNode(name=name,
                                                destination_compartments=[dest_symp, dest_asym],
                                                probability_params=params.probSym)

    # change nodes for treatment outcomes
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            for a in range(len(AB)):
                i = covert_symp_susp_antibio.get_row_index(symp_state=s, susp_profile=p, antibiotic=a)
                name = 'Tx outcome for ' + covert_symp_susp_antibio.get_str_symp_susp_antibio(
                    symp_state=s, susp_profile=p, antibiotic=a)

                # if drug-susceptible
                if p == RestProfile.SUS.value:
                    dest_succ = S
                    # failure due to the emergence of resistance
                    prob_fail = params.probResEmerge[a]
                    if a == AB.PEN.value:
                        dest_fail = ifs_symp_from_emerg_res[covert_symp_susp.get_row_index(
                            symp_state=s, susp_profile=RestProfile.PEN.value)]
                    elif a == AB.CFX.value:
                        dest_fail = ifs_symp_from_emerg_res[
                            covert_symp_susp.get_row_index(symp_state=s, susp_profile=RestProfile.CFX.value)]
                    else:
                        raise ValueError('Invalid antibiotic.')

                # if resistance to PEN
                elif p == RestProfile.PEN.value:
                    if a == AB.PEN.value:
                        # failure due to ineffective treatment
                        j = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
                        dest_succ = ifs_re_tx[j]
                        dest_fail = ifs_re_tx[j]
                        prob_fail = Constant(1)
                    elif a == AB.CFX:
                        # success or resistance
                        dest_succ = S
                        dest_fail = ifs_re_tx[covert_symp_susp.get_row_index(
                            symp_state=s, susp_profile=RestProfile.PEN_CFX.value)]
                        prob_fail = params.probResEmerge
                    else:
                        raise ValueError('Invalid antibiotic.')

                # if resistance to CFX
                elif p == RestProfile.CFX:
                    if a == AB.PEN:
                        # success or resistance
                        dest_succ = S
                        dest_fail = ifs_re_tx[covert_symp_susp.get_row_index(
                            symp_state=s, susp_profile=RestProfile.PEN_CFX.value)]
                        prob_fail = params.probResEmerge
                    elif a == AB.CFX:
                        # failure
                        j = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
                        dest_succ = ifs_re_tx[j]
                        dest_fail = ifs_re_tx[j]
                        prob_fail = Constant(1)
                    else:
                        raise ValueError('Invalid antibiotic.')

                elif p == RestProfile.PEN_CFX:
                    # failure
                    j = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
                    dest_succ = ifs_re_tx[j]
                    dest_fail = ifs_re_tx[j]
                    prob_fail = Constant(1)
                else:
                    raise ValueError('Invalid treatment outcome.')

                ifs_tx_outcomes[i] = ChanceNode(name=name,
                                                destination_compartments=[dest_fail, dest_succ],
                                                probability_params=prob_fail)

    # rapid tests
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            name = 'Rapid test in '+covert_symp_susp.get_str_symp_susp(symp_state=s, susp_profile=p)
            j_pos = covert_symp_susp_antibio.get_row_index(symp_state=s, susp_profile=p, antibiotic=AB.CFX.value)
            j_neg = covert_symp_susp_antibio.get_row_index(symp_state=s, susp_profile=p, antibiotic=AB.PEN.value)

            dest_pos = ifs_tx_outcomes[j_pos]
            dest_neg = ifs_tx_outcomes[j_neg]

            if s in (RestProfile.PEN.value or RestProfile.PEN_CFX.value):
                prob_pos = params.sens
            else:
                prob_pos = params.oneMinusSpec

            i = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
            ifs_rapid_tests[i] = ChanceNode(name=name,
                                            destination_compartments=[dest_pos, dest_neg],
                                            probability_params=prob_pos)

            if_rest_to[p].append(ifs_rapid_tests[i])

    # ------------- compartment histories ---------------
    # set up prevalence, incidence, and cumulative incidence to collect
    if sets.ifCollectTrajsOfCompartments:
        S.setup_history(collect_prev=True)
        for i in Is:
            i.setup_history(collect_prev=True)
        for f in Fs:
            f.setup_history(collect_prev=True)

    # ------------- summation statistics ---------------
    # population size
    all_comparts = [S]
    all_comparts.extend(Is)
    all_comparts.extend(Fs)

    pop_size = SumPrevalence(name='Population size', compartments=all_comparts)

    # number infected
    infected_comparts = [Is]
    infected_comparts.extend(Fs)
    n_infected = SumPrevalence(name='Number infected',
                               compartments=infected_comparts)
    prevalence = RatioTimeSeries(name='Prevalence',
                                 numerator_sum_time_series=n_infected,
                                 denominator_sum_time_series=pop_size,
                                 if_surveyed=True)

    # rate of gonorrhea cases
    new_cases = SumIncidence(name='Number of gonorrhea cases',
                             compartments=ifs_symp_from_S)
    gono_rate = RatioTimeSeries(name='Rate of gonorrhea cases',
                                numerator_sum_time_series=new_cases,
                                denominator_sum_time_series=pop_size,
                                if_surveyed=True)

    # % cases symptomatic
    cases_symptomatic = SumIncidence(name='Number of cases symptomatic',
                                     compartments=Is_symp)
    percent_cases_symptomatic = RatioTimeSeries(name='Proportion of cases symptomatic',
                                                numerator_sum_time_series=cases_symptomatic,
                                                denominator_sum_time_series=new_cases,
                                                if_surveyed=True)

    # cases resistant to drug A
    cases_resistant_by_profile = []
    perc_cases_resistant_by_profile = []

    for p in range(len(RestProfile)):
        cases = SumIncidence(name='Number of cases resistant to ' + REST_PROFILES[p],
                             compartments=if_rest_to[p])
        perc_cases = RatioTimeSeries(
            name='Proportion of cases resistant to ' + REST_PROFILES[p],
            numerator_sum_time_series=cases,
            denominator_sum_time_series=new_cases,
            if_surveyed=True,
            survey_size_param=params.surveySize)

        cases_resistant_by_profile.append(cases)
        perc_cases_resistant_by_profile.append(perc_cases)

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

    if False:
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
    inf_events = []
    for p in range(len(RestProfile)):
        S.add_event(EpiDepEvent(
                name='Infection | ' + REST_PROFILES[p],
                destination=ifs_symp_from_S[p],
                generating_pathogen=p))

    # add events to infection compartments
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            i = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
            compart_name = Is[i].name

            # natural recovery
            Is[i].add_event(EpiIndepEvent(
                name='Natural recovery | ' + compart_name,
                rate_param=params.rateNaturalRecovery,
                destination=S))
            # screened
            Is[i].add_event(EpiIndepEvent(
                name='Screening | ' + compart_name,
                rate_param=params.rateScreened,
                destination=ifs_rapid_tests[i]))
            # seeking treatment
            if s == SympStat.SYMP.value:
                Is[i].add_event(EpiIndepEvent(
                    name='Seeking treatment | ' + compart_name,
                    rate_param=params.rateTreatment,
                    destination=ifs_rapid_tests[i]))

    # add events to infection compartments after treatment failure
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            i = covert_symp_susp.get_row_index(symp_state=s, susp_profile=p)
            compart_name = Is[i].name
            # treatment
            Fs[i].add_event(EpiIndepEvent(
                name='Re-Tx | ' + compart_name,
                rate_param=params.rateTreatment,
                destination=S
            ))

    # ------------- population the model ---------------
    # populate the model
    chance_nodes = ifs_symp_from_S
    chance_nodes.extend(ifs_re_tx)
    chance_nodes.extend(ifs_symp_from_emerg_res)
    chance_nodes.extend(ifs_tx_outcomes)
    chance_nodes.extend(ifs_rapid_tests)

    list_of_sum_time_series = [pop_size, n_infected, new_cases, cases_symptomatic]
    list_of_sum_time_series.extend(cases_resistant_by_profile)

    list_of_ratio_time_series = [prevalence, gono_rate, percent_cases_symptomatic]
    list_of_ratio_time_series.extend(perc_cases_resistant_by_profile)

    model.populate(compartments=all_comparts,
                   chance_nodes=chance_nodes,
                   list_of_sum_time_series=list_of_sum_time_series,
                   list_of_ratio_time_series=list_of_ratio_time_series,
                   interventions=None, # [first_line_tx_with_A, first_line_tx_with_M],
                   features=None, # [f_perc_resist_A, f_if_M_ever_switched_on],
                   conditions=None, # [A_in_condition, A_out_condition, M_is_never_used, ],
                   parameters=params)
