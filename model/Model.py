from SimPy.Parameters import Constant
from apace.CalibrationSupport import FeasibleConditions
from apace.Control import InterventionAffectingEvents, ConditionBasedDecisionRule
from apace.FeaturesAndConditions import FeatureSurveillance, FeatureIntervention, \
    ConditionOnFeatures, ConditionOnConditions, ConditionAlwaysFalse
from apace.ModelObjects import Compartment, ChanceNode, EpiIndepEvent, EpiDepEvent
from apace.TimeSeries import SumPrevalence, SumIncidence, RatioTimeSeries
from definitions import RestProfile, AB, SympStat, REST_PROFILES, TreatmentOutcome, \
    ConvertSympAndSuspAndAntiBio, get_profile_after_resit_or_failure
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
    ifs_rapid_tests = [None] * covert_symp_susp.length
    ifs_re_tx = [None] * covert_symp_susp.length
    ifs_symp_from_emerg_rest = [None] * len(RestProfile)
    ifs_tx_outcomes = [None] * covert_symp_susp_antibio.length
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
    # counting M used
    counting_tx_M = ChanceNode(name='Tx with M',
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
            # profile after the rise of resistance to CFX
            next_p, reason_for_failure = get_profile_after_resit_or_failure(
                rest_profile=p, antibiotic=AB.CFX)

            # if treatment failure is because of the development of resistance
            if reason_for_failure == TreatmentOutcome.RESISTANCE:

                # name of this compartment
                str_symp_susp = covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p)
                i = covert_symp_susp.get_row_index(symp_state=s, rest_profile=p)

                # destinations
                dest_rest = Fs[covert_symp_susp.get_row_index(symp_state=s, rest_profile=next_p)]
                dest_succ = counting_success_CIP_TET_CFX
                ifs_resist_after_re_tx_cfx[i] = ChanceNode(
                    name='If resistance after re-Tx-CFX in '+str_symp_susp,
                    destination_compartments=[dest_rest, dest_succ],
                    probability_params=params.probResEmerge[AB.CFX.value])

    # if re-tx after ineffective treatment or resistance development
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'If re-tx after ineffective tx/resistance development in ' \
                   + covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_susp.get_row_index(symp_state=s, rest_profile=p)

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

    # chance nodes for if symptomatic after emergence of resistance during treatment
    for p in range(n_rest_profiles):
        name = 'If symptomatic after moving to {} due to the emergence of resistance'.format(
            covert_symp_susp.get_str_susp(susp_profile=p))

        dest_symp = ifs_re_tx[covert_symp_susp.get_row_index(symp_state=SympStat.SYMP.value, rest_profile=p)]
        dest_asym = ifs_re_tx[covert_symp_susp.get_row_index(symp_state=SympStat.ASYM.value, rest_profile=p)]
        ifs_symp_from_emerg_rest[p] = ChanceNode(name=name,
                                                 destination_compartments=[dest_symp, dest_asym],
                                                 probability_params=params.probSym)

    # rapid test (for now everyone will receive CFX)
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'Rapid test (all receive CFX) in ' \
                   + covert_symp_susp.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_susp.get_row_index(symp_state=s, rest_profile=p)

            # treatment outcomes
            next_p, reason_for_failure = get_profile_after_resit_or_failure(
                rest_profile=p, antibiotic=AB.CFX)

            # if treatment failure is because of the development of resistance
            if reason_for_failure == TreatmentOutcome.RESISTANCE:
                # destinations
                dest_rest = ifs_symp_from_emerg_rest[
                    covert_symp_susp.get_row_index(symp_state=s, rest_profile=next_p)]
                dest_succ = counting_success_CIP_TET_CFX
                ifs_rapid_tests[i] = ChanceNode(
                    name=name,
                    destination_compartments=[dest_rest, dest_succ],
                    probability_params=params.probResEmerge[AB.CFX.value])

            elif reason_for_failure == TreatmentOutcome.INEFFECTIVE:
                ifs_rapid_tests[i] = ChanceNode(
                    name=name,
                    destination_compartments=[Fs[i], Fs[i]],
                    probability_params=Constant(1))

    # ------------- compartment histories ---------------
    # set up prevalence, incidence, and cumulative incidence to collect
    if sets.ifCollectTrajsOfCompartments:
        S.setup_history(collect_prev=True)
        for i in Is:
            i.setup_history(collect_prev=True)
        for f in Fs:
            f.setup_history(collect_prev=True)
        # for r in ifs_CIP_rapid_tests:
        #     r.setup_history(collect_incd=True)
        for t in ifs_tx_outcomes:
            t.setup_history(collect_incd=True)
    counting_success_CIP_TET_CFX.setup_history(collect_incd=True)
    counting_tx_M.setup_history(collect_incd=True)

    # ------------- summation statistics ---------------
    # population size
    all_comparts = [S] + Is + Fs
    pop_size = SumPrevalence(name='Population size', compartments=all_comparts)

    # number infected
    n_infected = SumPrevalence(name='Infected',
                               compartments=Is + Fs)
    prevalence = RatioTimeSeries(name='Prevalence',
                                 numerator_sum_time_series=n_infected,
                                 denominator_sum_time_series=pop_size,
                                 if_surveyed=True)

    # rate of new gonorrhea cases
    n_cases = SumIncidence(name='New cases',
                           compartments=ifs_CIP_rapid_tests)
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
        compartments=[counting_success_CIP_TET_CFX, counting_tx_M] + ifs_counting_tx_M)
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
                destination=ifs_CIP_rapid_tests[i],
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
                    destination=ifs_CIP_rapid_tests[i],
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
                dest = counting_tx_M
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
    chance_nodes.extend(ifs_CIP_rapid_tests)
    chance_nodes.extend(ifs_counting_tx_M)
    chance_nodes.extend([counting_success_CIP_TET_CFX, counting_tx_M])

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
