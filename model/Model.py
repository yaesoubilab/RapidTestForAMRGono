from SimPy.Parameters import Constant
from apace.CalibrationSupport import FeasibleConditions
from apace.Control import InterventionAffectingEvents, ConditionBasedDecisionRule
from apace.FeaturesAndConditions import FeatureSurveillance, FeatureIntervention, \
    ConditionOnFeatures, ConditionOnConditions, ConditionAlwaysFalse
from apace.ModelObjects import Compartment, ChanceNode, EpiIndepEvent, EpiDepEvent
from apace.TimeSeries import SumPrevalence, SumIncidence, RatioTimeSeries
from definitions import RestProfile, AB, SympStat, REST_PROFILES, ConvertSympAndSuspAndAntiBio
from model.ModelParameters import Parameters


def build_model(model):

    # constants
    n_symp_states = len(SympStat)
    n_rest_profiles = len(RestProfile)

    # model settings
    sets = model.settings

    # object to convert (resistance profile, symptom status) to row index or string
    covert_symp_susp = ConvertSympAndSuspAndAntiBio(
        n_symp_stats=len(SympStat), n_rest_profiles=len(RestProfile))
    # object to convert (resistance profile, symptom status, antibiotic) to row index or string
    covert_symp_susp_antibio = ConvertSympAndSuspAndAntiBio(
        n_symp_stats=len(SympStat), n_rest_profiles=len(RestProfile), n_antibiotics=len(AB))

    # model parameters
    params = Parameters(model_sets=sets)

    # ------------- model compartments ---------------
    Is = [None] * covert_symp_susp.length
    Fs = [None] * covert_symp_susp.length
    ifs_symp_from_S = [None] * len(RestProfile)
    ifs_re_tx = [None] * covert_symp_susp.length
    ifs_symp_from_emerg_rest = [None] * len(RestProfile)
    ifs_tx_outcomes = [None] * covert_symp_susp_antibio.length
    ifs_rapid_tests = [None] * covert_symp_susp.length
    ifs_counting_tx_M = [None] * covert_symp_susp.length
    ifs_resist_after_re_tx_cfx = [None] * covert_symp_susp.length

    ifs_symp = []   # chance nodes counting symptomatic cases
    ifs_rest_to = []  # chance nodes counting resistance or intermediate resistance to different profiles
    for p in RestProfile:
        ifs_rest_to.append([])

    # S
    S = Compartment(name='S', size_par=params.sizeS,
                    susceptibility_params=[Constant(value=1)]*len(RestProfile))

    # Is (infectious) and Fs (infectious after treatment failure)
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            # name of this compartment
            str_symp_susp = covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_susp.get_row_index(symp_state=s, rest_profile=p)

            # infectivity parameters
            infectivity_params = [Constant(value=0)] * covert_symp_susp.nRestProfiles
            infectivity_params[p] = params.infectivityByRestProfile[p]

            # Is
            Is[i] = Compartment(name='I ' + str_symp_susp,
                                size_par=params.sizeIBySympAndRest[i],
                                if_empty_to_eradicate=True,
                                infectivity_params=infectivity_params)
            # Fs: infectious compartments after treatment failure
            Fs[i] = Compartment(name='F ' + str_symp_susp,
                                if_empty_to_eradicate=True,
                                infectivity_params=infectivity_params)

    # counting successful treatments with any drug
    counting_success_CIP_TET_CFX = ChanceNode(
        name='Successful Tx with CIP, TET, or CFX',
        destination_compartments=[S, S],
        probability_params=Constant(1))
    # counting successful re-treatments with CFX
    counting_success_re_tx_M = ChanceNode(name='Successful re-Tx with M',
                                          destination_compartments=[S, S],
                                          probability_params=Constant(1))

    # if symptomatic after infection
    for p in range(n_rest_profiles):
        name = 'If symptomatic to ' + covert_symp_susp.get_str_susp(susp_profile=p)
        dest_symp = Is[covert_symp_susp.get_row_index(symp_state=SympStat.SYMP.value, rest_profile=p)]
        dest_asym = Is[covert_symp_susp.get_row_index(symp_state=SympStat.ASYM.value, rest_profile=p)]
        ifs_symp_from_S[p] = ChanceNode(name=name,
                                        destination_compartments=[dest_symp, dest_asym],
                                        probability_params=params.probSym)

    # if resistance emerges after re-tx with CFX
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            # profiles for which resistance to CFX might arise
            if p in (RestProfile.CIP.value, RestProfile.TET.value, RestProfile.CIP_TET):

                # name of this compartment
                str_symp_susp = covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p)
                i = covert_symp_susp.get_row_index(symp_state=s, rest_profile=p)

                # profile after the rise of resistance to CFX
                if p == RestProfile.CIP.value:
                    next_p = RestProfile.CIP_CFX
                elif p == RestProfile.TET.value:
                    next_p = RestProfile.TET_CFX
                elif p == RestProfile.CIP_TET:
                    next_p = RestProfile.CIP_TET_CFX
                else:
                    raise Exception('Invalid value of resitance profile.')

                dest_rest = Fs[covert_symp_susp.get_row_index(symp_state=s, rest_profile=next_p)]
                dest_succ = counting_success_CIP_TET_CFX
                ifs_resist_after_re_tx_cfx[i] = ChanceNode(
                    name='If resistance after re-Tx-CFX in '+str_symp_susp,
                    destination_compartments=[dest_rest, dest_succ],
                    probability_params=params.probResEmerge[AB.CFX.value])

    # if re-tx after ineffective treatment or resistance development
    i = 0
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'If re-tx after ineffective tx/resistance development in ' \
                   + covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p)
            dest_yes = Fs[i]
            dest_no = Is[i]
            # symptomatic cases will always seek re-treatment and asymptomatic cases never
            if s == SympStat.SYMP.value:
                prob_yes = Constant(1)
            elif s == SympStat.ASYM.value:
                prob_yes = Constant(0)
            else:
                raise ValueError('Invalid symptom status.')

            ifs_re_tx[i] = ChanceNode(name=name,
                                      destination_compartments=[dest_yes, dest_no],
                                      probability_params=prob_yes)
            i += 1

    # chance nodes for if symptomatic after emergence of resistance during treatment
    for p in range(n_rest_profiles):
        name = 'If symptomatic after moving to {} due to the emergence of resistance'.format(
            covert_symp_susp.get_str_susp(susp_profile=p))
        dest_symp = ifs_re_tx[covert_symp_susp.get_row_index(symp_state=SympStat.SYMP.value, rest_profile=p)]
        dest_asym = ifs_re_tx[covert_symp_susp.get_row_index(symp_state=SympStat.ASYM.value, rest_profile=p)]
        ifs_symp_from_emerg_rest[p] = ChanceNode(name=name,
                                                 destination_compartments=[dest_symp, dest_asym],
                                                 probability_params=params.probSym)

    # chance nodes for treatment outcomes
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            for a in range(len(AB)):
                i = covert_symp_susp_antibio.get_row_index(symp_state=s, rest_profile=p, antibiotic=a)
                name = 'Tx outcome for ' + covert_symp_susp_antibio.get_str_symp_rest_antibio(
                    symp_state=s, rest_profile=p, antibiotic=a)

                # if drug-susceptible
                if p == RestProfile.SUS.value:
                    # failure due to the emergence of resistance
                    dest_succ = counting_success_CIP_TET_CFX
                    prob_fail = params.probResEmerge[a]

                    # if resistance emerges while receiving PEN
                    if a == AB.PEN.value:
                        # decide where to go next depending on the symptom status
                        if s == SympStat.SYMP.value:
                            dest_fail = ifs_re_tx[
                                covert_symp_susp.get_row_index(symp_state=s, rest_profile=RestProfile.PEN.value)]
                        elif s == SympStat.ASYM.value:
                            dest_fail = ifs_symp_from_emerg_rest[RestProfile.PEN.value]
                        else:
                            raise ValueError('Invalid symptom status.')

                    # if resistance emerges while receiving CFX
                    elif a == AB.CFX.value:
                        # decide where to go next depending on the symptom status
                        if s == SympStat.SYMP.value:
                            dest_fail = ifs_re_tx[
                                covert_symp_susp.get_row_index(symp_state=s, rest_profile=RestProfile.CFX.value)]
                        elif s == SympStat.ASYM.value:
                            dest_fail = ifs_symp_from_emerg_rest[RestProfile.CFX.value]
                        else:
                            raise ValueError('Invalid symptom status.')
                    else:
                        raise ValueError('Invalid antibiotic.')

                # if resistance to PEN or CFX
                elif p in (RestProfile.PEN.value, RestProfile.CFX.value):
                    # failure due to ineffective treatment
                    if (p == RestProfile.PEN.value and a == AB.PEN.value) or \
                            (p == RestProfile.CFX.value and a == AB.CFX.value):
                        # failure due to ineffective treatment
                        j = covert_symp_susp.get_row_index(symp_state=s, rest_profile=p)
                        dest_succ = ifs_re_tx[j]
                        dest_fail = ifs_re_tx[j]
                        prob_fail = Constant(1)

                    # failure due to the emergence of resistance
                    elif (p == RestProfile.PEN.value and a == AB.CFX.value) or \
                            (p == RestProfile.CFX.value and a == AB.PEN.value):
                        # success or resistance
                        dest_succ = counting_success_CIP_TET_CFX
                        prob_fail = params.probResEmerge[a]

                        # decide where to go next depending on the symptom status
                        if s == SympStat.SYMP.value:
                            dest_fail = ifs_re_tx[
                                covert_symp_susp.get_row_index(symp_state=s, rest_profile=RestProfile.PEN_CFX.value)]
                        elif s == SympStat.ASYM.value:
                            dest_fail = ifs_symp_from_emerg_rest[RestProfile.PEN_CFX.value]
                        else:
                            raise ValueError('Invalid symptom status.')
                    else:
                        raise ValueError('Invalid antibiotic.')

                # if resistance to PEN and CFX
                elif p == RestProfile.PEN_CFX.value:
                    # failure
                    j = covert_symp_susp.get_row_index(symp_state=s, rest_profile=p)
                    dest_succ = ifs_re_tx[j]
                    dest_fail = ifs_re_tx[j]
                    prob_fail = Constant(1)
                else:
                    raise ValueError('Invalid treatment outcome.')

                ifs_tx_outcomes[i] = ChanceNode(name=name,
                                                destination_compartments=[dest_fail, dest_succ],
                                                probability_params=prob_fail)

    # rapid tests
    i = 0
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):

            # indices of destination
            # if results is positive (i.e., susceptibility to PEN)
            j_pos = covert_symp_susp_antibio.get_row_index(symp_state=s, rest_profile=p, antibiotic=AB.PEN.value)
            # if results is negative (i.e., non-susceptibility to PEN)
            j_neg = covert_symp_susp_antibio.get_row_index(symp_state=s, rest_profile=p, antibiotic=AB.CFX.value)

            # destinations
            dest_pos = ifs_tx_outcomes[j_pos]
            dest_neg = ifs_tx_outcomes[j_neg]

            # probability that the test result is positive
            if p in (RestProfile.PEN.value, RestProfile.PEN_CFX.value):
                prob_pos = params.oneMinusSpec
            else:
                prob_pos = params.sens

            ifs_rapid_tests[i] = ChanceNode(name='Rapid test in ' +
                                                 covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p),
                                            destination_compartments=[dest_pos, dest_neg],
                                            probability_params=prob_pos)
            ifs_counting_tx_M[i] = ChanceNode(name='Counting Tx-M from ' +
                                                   covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p),
                                              destination_compartments=[S, S],
                                              probability_params=Constant(1))
            # if symptomatic
            if s == SympStat.SYMP.value:
                ifs_symp.append(ifs_rapid_tests[i])
                ifs_symp.append(ifs_counting_tx_M[i])
            # by resistance profile
            ifs_rest_to[p].append(ifs_rapid_tests[i])
            ifs_rest_to[p].append(ifs_counting_tx_M[i])
            i += 1

    # ------------- compartment histories ---------------
    # set up prevalence, incidence, and cumulative incidence to collect
    if sets.ifCollectTrajsOfCompartments:
        S.setup_history(collect_prev=True)
        for i in Is:
            i.setup_history(collect_prev=True)
        for f in Fs:
            f.setup_history(collect_prev=True)
        for r in ifs_rapid_tests:
            r.setup_history(collect_incd=True)
        for r in ifs_counting_tx_M:
            r.setup_history(collect_incd=True)
        for t in ifs_tx_outcomes:
            t.setup_history(collect_incd=True)
        counting_success_CIP_TET_CFX.setup_history(collect_incd=True)
        counting_success_re_tx_M.setup_history(collect_incd=True)

    # ------------- summation statistics ---------------
    # population size
    all_comparts = [S]
    all_comparts.extend(Is)
    all_comparts.extend(Fs)

    pop_size = SumPrevalence(name='Population size', compartments=all_comparts)

    # number infected
    infected_comparts = Is
    infected_comparts.extend(Fs)
    n_infected = SumPrevalence(name='Infected',
                               compartments=infected_comparts)
    prevalence = RatioTimeSeries(name='Prevalence',
                                 numerator_sum_time_series=n_infected,
                                 denominator_sum_time_series=pop_size,
                                 if_surveyed=True)

    # rate of new gonorrhea cases
    n_cases = SumIncidence(name='New cases',
                           compartments=ifs_rapid_tests + ifs_counting_tx_M)
    gono_rate = RatioTimeSeries(name='Rate of gonorrhea cases',
                                numerator_sum_time_series=n_cases,
                                denominator_sum_time_series=pop_size,
                                if_surveyed=True,
                                collect_stat_after_warm_up=True)

    # % cases symptomatic
    n_cases_sympt = SumIncidence(name='New cases symptomatic',
                                 compartments=ifs_symp)
    perc_cases_sympt = RatioTimeSeries(name='Proportion of cases symptomatic',
                                       numerator_sum_time_series=n_cases_sympt,
                                       denominator_sum_time_series=n_cases,
                                       if_surveyed=True)

    # cases by resistance profile
    n_cases_by_resistance_profile = []
    perc_cases_by_resistance_profile = []
    for p in range(n_rest_profiles):
        n_resistant_cases = SumIncidence(name='Cases ' + REST_PROFILES[p],
                                         compartments=ifs_rest_to[p])
        perc_cases_resistant = RatioTimeSeries(
            name='Proportion of cases resistant to ' + REST_PROFILES[p],
            numerator_sum_time_series=n_resistant_cases,
            denominator_sum_time_series=n_cases,
            if_surveyed=True,
            survey_size_param=params.surveySize)

        n_cases_by_resistance_profile.append(n_resistant_cases)
        perc_cases_by_resistance_profile.append(perc_cases_resistant)

    # cases by resistance to specific CFX (including resistance to both PEN and CFX)
    n_cases_CFX_R = SumIncidence(
        name='Cases CFX-R or CFX+PEN-R',
        compartments=ifs_rest_to[RestProfile.CFX.value] + ifs_rest_to[RestProfile.PEN_CFX.value])
    perc_cases_CFX_R = RatioTimeSeries(
            name='Proportion of cases CFX-R or CFX+PEN-R',
            numerator_sum_time_series=n_cases_CFX_R,
            denominator_sum_time_series=n_cases,
            if_surveyed=True,
            survey_size_param=params.surveySize)

    # treated with any antibiotics
    n_treated = SumIncidence(
        name='Cases treated',
        compartments=[counting_success_CIP_TET_CFX, counting_success_re_tx_M] + ifs_counting_tx_M)
    n_treated_PEN_or_CFX = SumIncidence(name='Treated with PEN or CFX',
                                        compartments=[counting_success_CIP_TET_CFX])
    perc_treated_with_PEN_or_CFX = RatioTimeSeries(
        name='Proportion of cases treated with PEN or CFX',
        numerator_sum_time_series=n_treated_PEN_or_CFX,
        denominator_sum_time_series=n_treated,
        collect_stat_after_warm_up=True)

    # ------------- calibration targets ---------------
    if sets.calcLikelihood:
        # prevalence
        prevalence.add_feasible_conditions(feasible_conditions=FeasibleConditions(
            feasible_min=0.025, feasible_max=0.065))
        prevalence.add_calibration_targets(ratios=sets.prevMean,
                                           survey_sizes=sets.prevN)
        # gonorrhea rate
        gono_rate.add_feasible_conditions(feasible_conditions=FeasibleConditions(
            feasible_min=0.045, feasible_max=0.085))
        gono_rate.add_calibration_targets(ratios=sets.gonoRateMean,
                                          survey_sizes=sets.gonoRateN)
        # % cases symptomatic
        perc_cases_sympt.add_feasible_conditions(feasible_conditions=FeasibleConditions(
            feasible_min=0.5, feasible_max=1))
        perc_cases_sympt.add_calibration_targets(ratios=sets.percSympMean,
                                                 survey_sizes=sets.percSympN)

    # ------------- interventions ---------------
    # interventions
    first_line_tx_with_CFX = InterventionAffectingEvents(name='1st line therapy with CFX')
    first_line_tx_with_M = InterventionAffectingEvents(name='1st line therapy with M')

    # ------------- features ---------------
    # features
    f_perc_resist_CFX = FeatureSurveillance(name='Surveyed % of cases resistant to CFX',
                                            ratio_time_series_with_surveillance=perc_cases_CFX_R)
    f_if_M_ever_switched_on = FeatureIntervention(name='If Drug M ever switched on',
                                                  intervention=first_line_tx_with_M,
                                                  feature_type='if ever switched on')

    # ------------- conditions ---------------
    # conditions
    threshold = 0.05
    CFX_in_condition = ConditionOnFeatures(name='If % resistant to CFX is below threshold',
                                           features=[f_perc_resist_CFX],
                                           signs=['l'],
                                           thresholds=[threshold])
    CFX_out_condition = ConditionOnFeatures(name='If % resistant to CFX passes threshold',
                                            features=[f_perc_resist_CFX],
                                            signs=['ge'],
                                            thresholds=[threshold])
    M_is_never_used = ConditionOnFeatures(name='If M is ever used as 1st-line therapy',
                                          features=[f_if_M_ever_switched_on],
                                          signs=['e'],
                                          thresholds=[0])

    turn_on_tx_with_CFX = ConditionOnConditions(name='Turn on Tx-CFX',
                                                conditions=[CFX_in_condition, M_is_never_used])

    # ------------- decision rules ---------------
    # add decision rules to interventions
    first_line_tx_with_CFX.add_decision_rule(
        decision_rule=ConditionBasedDecisionRule(default_switch_value=1,
                                                 condition_to_turn_on=turn_on_tx_with_CFX,
                                                 condition_to_turn_off=CFX_out_condition))
    first_line_tx_with_M.add_decision_rule(
        decision_rule=ConditionBasedDecisionRule(default_switch_value=0,
                                                 condition_to_turn_on=CFX_out_condition,
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
    i = 0
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            compart_name = Is[i].name

            # natural recovery
            Is[i].add_event(EpiIndepEvent(
                name='Natural recovery | ' + compart_name,
                rate_param=params.rateNaturalRecovery,
                destination=S))
            # screened
            Is[i].add_event(EpiIndepEvent(
                name='Screening then CFX| ' + compart_name,
                rate_param=params.rateScreened,
                destination=ifs_rapid_tests[i],
                interv_to_activate=first_line_tx_with_CFX))
            Is[i].add_event(EpiIndepEvent(
                name='Screening then M| ' + compart_name,
                rate_param=params.rateScreened,
                destination=ifs_counting_tx_M[i],
                interv_to_activate=first_line_tx_with_M))
            # seeking treatment
            if s == SympStat.SYMP.value:
                Is[i].add_event(EpiIndepEvent(
                    name='Seeking treatment then CFX | ' + compart_name,
                    rate_param=params.rateTreatment,
                    destination=ifs_rapid_tests[i],
                    interv_to_activate=first_line_tx_with_CFX))
                Is[i].add_event(EpiIndepEvent(
                    name='Seeking treatment then M | ' + compart_name,
                    rate_param=params.rateTreatment,
                    destination=ifs_counting_tx_M[i],
                    interv_to_activate=first_line_tx_with_M))
            i += 1

    # add events to infection compartments after treatment failure
    i = 0
    for s in range(len(SympStat)):
        for p in range(len(RestProfile)):
            compart_name = Fs[i].name
            # destination
            if p == RestProfile.PEN.value:
                dest = ifs_resist_after_re_tx_cfx[s]
            else:
                dest = counting_success_re_tx_M
            # treatment
            Fs[i].add_event(EpiIndepEvent(
                name='Re-Tx | ' + compart_name,
                rate_param=params.rateRetreatment,
                destination=dest
            ))
            i += 1

    # ------------- populate the model ---------------
    # populate the model
    chance_nodes = ifs_symp_from_S
    chance_nodes.extend(ifs_re_tx)
    chance_nodes.extend(ifs_symp_from_emerg_rest)
    chance_nodes.extend(ifs_tx_outcomes)
    chance_nodes.extend(ifs_rapid_tests)
    chance_nodes.extend(ifs_counting_tx_M)
    chance_nodes.extend([counting_success_CIP_TET_CFX, counting_success_re_tx_M])

    list_of_sum_time_series = [pop_size, n_infected, n_cases, n_cases_CFX_R,
                               n_cases_sympt, n_treated, n_treated_PEN_or_CFX]
    list_of_sum_time_series.extend(n_cases_by_resistance_profile)

    list_of_ratio_time_series = [prevalence, gono_rate, perc_cases_CFX_R,
                                 perc_cases_sympt, perc_treated_with_PEN_or_CFX]
    list_of_ratio_time_series.extend(perc_cases_by_resistance_profile)

    model.populate(compartments=all_comparts,
                   chance_nodes=chance_nodes,
                   list_of_sum_time_series=list_of_sum_time_series,
                   list_of_ratio_time_series=list_of_ratio_time_series,
                   interventions=[first_line_tx_with_CFX, first_line_tx_with_M],
                   features=[f_perc_resist_CFX, f_if_M_ever_switched_on],
                   conditions=[CFX_in_condition, CFX_out_condition, M_is_never_used, turn_on_tx_with_CFX],
                   parameters=params)
