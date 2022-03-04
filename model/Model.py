from SimPy.Parameters import Constant
from apace.CalibrationSupport import FeasibleConditions
from apace.Control import InterventionAffectingEvents
from apace.ModelObjects import Compartment, ChanceNode, EpiIndepEvent, EpiDepEvent
from apace.TimeSeries import SumPrevalence, SumIncidence, RatioTimeSeries
from definitions import RestProfile, AB, SympStat, REST_PROFILES, ANTIBIOTICS, TreatmentOutcome, \
    ConvertSympAndResitAndAntiBio, get_profile_after_resit_or_failure
from model.ModelParameters import Parameters


def build_model(model):

    # constants
    n_symp_states = len(SympStat)
    n_rest_profiles = len(RestProfile)

    # model settings
    sets = model.settings

    # object to convert (resistance profile, symptom status) to row index or string
    covert_symp_resist = ConvertSympAndResitAndAntiBio(
        n_symp_stats=len(SympStat), n_rest_profiles=len(RestProfile))
    # object to convert (resistance profile, symptom status, antibiotic) to row index or string
    covert_symp_resist_antibio = ConvertSympAndResitAndAntiBio(
        n_symp_stats=len(SympStat), n_rest_profiles=len(RestProfile), n_antibiotics=len(AB))

    # model parameters
    params = Parameters(model_sets=sets)

    # ------------- model compartments ---------------
    Is = [None] * covert_symp_resist.length
    Fs = [None] * covert_symp_resist.length
    ifs_symp_from_S = [None] * len(RestProfile)

    # chance nodes to model rapid susceptibility test for CIP and TET
    ifs_will_receive_rapid_test = [None] * covert_symp_resist.length
    ifs_rapid_CIP_outcome = [None] * covert_symp_resist.length
    ifs_rapid_TET_outcome_after_susp_CIP = [None] * covert_symp_resist.length
    ifs_rapid_TET_outcome_after_reduced_susp_CIP = [None] * covert_symp_resist.length

    # chance nodes to model if using CIP or TET when susceptible to both
    ifs_CIP_or_TET = [None] * covert_symp_resist.length
    # chance nodes to model treatment outcome by resistance profile, symptom status, and antibiotic for tx
    ifs_tx_outcome = [None] * covert_symp_resist_antibio.length

    # chance nodes to model if will seek re-treatment by resistance profile and symptom status
    ifs_re_tx = [None] * covert_symp_resist.length
    ifs_symp_from_emerg_rest = [None] * len(RestProfile)
    ifs_resist_after_re_tx_cfx = [None] * covert_symp_resist.length

    counting_symp = []   # chance nodes counting symptomatic cases
    counting_rest_to = []  # chance nodes counting resistance or intermediate resistance to different profiles
    for p in RestProfile:
        counting_rest_to.append([])

    # S
    S = Compartment(name='S', size_par=params.sizeS,
                    susceptibility_params=[Constant(value=1)]*len(RestProfile))

    # Is (infectious) and Fs (infectious after treatment failure)
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            # name of this compartment
            str_symp_susp = covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

            # infectivity parameters
            infectivity_params = [Constant(value=0)] * covert_symp_resist.nRestProfiles
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
    # counting successful Tx with each antibiotic
    counting_tx_success_by_ab = [None] * len(AB)
    for a in range(len(AB)):
        counting_tx_success_by_ab[a] = ChanceNode(
            name='Counting success tx with ' + ANTIBIOTICS[a],
            destination_compartments=[counting_success_CIP_TET_CFX, counting_success_CIP_TET_CFX],
            probability_params=Constant(1))

    # if symptomatic after infection
    for p in range(n_rest_profiles):
        name = 'If symptomatic to ' + covert_symp_resist.get_str_susp(susp_profile=p)
        dest_symp = Is[covert_symp_resist.get_row_index(symp_state=SympStat.SYMP.value, rest_profile=p)]
        dest_asym = Is[covert_symp_resist.get_row_index(symp_state=SympStat.ASYM.value, rest_profile=p)]
        ifs_symp_from_S[p] = ChanceNode(name=name,
                                        destination_compartments=[dest_symp, dest_asym],
                                        probability_params=params.probSym)

    # if resistance emerges after re-tx with CFX
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):

            # profile after the rise of resistance to CFX
            next_p, reason_for_failure = get_profile_after_resit_or_failure(
                rest_profile=p, antibiotic=AB.CFX)

            # name of this compartment
            str_symp_susp = covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

            # if treatment failure is because of the development of resistance
            if reason_for_failure == TreatmentOutcome.RESISTANCE:

                # destinations
                dest_rest = Fs[covert_symp_resist.get_row_index(symp_state=s, rest_profile=next_p)]
                dest_succ = counting_tx_success_by_ab[AB.CFX.value]
                ifs_resist_after_re_tx_cfx[i] = ChanceNode(
                    name='If resistance after re-tx with CFX in '+str_symp_susp,
                    destination_compartments=[dest_rest, dest_succ],
                    probability_params=params.probResEmerge[AB.CFX.value])
            else:
                # destinations
                dest_fail = counting_tx_M
                dest_succ = counting_tx_success_by_ab[AB.CFX.value]
                # this is a dummy chance node which will never be used.
                ifs_resist_after_re_tx_cfx[i] = ChanceNode(
                    name='If treatment failure after re-tx with CFX in ' + str_symp_susp,
                    destination_compartments=[dest_fail, dest_succ],
                    probability_params=Constant(1))

    # if re-tx after ineffective treatment or resistance development
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'If re-tx after ineff tx/rest dev in ' \
                   + covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

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

    # chance nodes for if symptomatic after emergence of resistance after 1st line treatment
    for p in range(n_rest_profiles):
        name = 'If symp after developing resistance to {}'.format(
            covert_symp_resist.get_str_susp(susp_profile=p))

        dest_symp = ifs_re_tx[covert_symp_resist.get_row_index(symp_state=SympStat.SYMP.value, rest_profile=p)]
        dest_asym = ifs_re_tx[covert_symp_resist.get_row_index(symp_state=SympStat.ASYM.value, rest_profile=p)]
        ifs_symp_from_emerg_rest[p] = ChanceNode(name=name,
                                                 destination_compartments=[dest_symp, dest_asym],
                                                 probability_params=params.probSym)

    # treatment outcomes
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            for a in range(len(AB)):

                name = 'Outcome of ' + covert_symp_resist_antibio.get_str_symp_rest_antibio(
                    symp_state=s, rest_profile=p, antibiotic=a)
                i = covert_symp_resist_antibio.get_row_index(symp_state=s, rest_profile=p, antibiotic=a)

                # treatment outcomes
                next_p, reason_for_failure = get_profile_after_resit_or_failure(
                    rest_profile=p, antibiotic=AB(a))

                next_i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=next_p)

                # if treatment failure is because of the development of resistance
                if reason_for_failure == TreatmentOutcome.RESISTANCE:

                    # destinations
                    if s == SympStat.SYMP.value:
                        dest_rest = ifs_re_tx[next_i]
                    else:
                        dest_rest = ifs_symp_from_emerg_rest[next_p]

                    dest_succ = counting_tx_success_by_ab[a]
                    ifs_tx_outcome[i] = ChanceNode(
                        name=name,
                        destination_compartments=[dest_rest, dest_succ],
                        probability_params=params.probResEmerge[a])

                elif reason_for_failure == TreatmentOutcome.INEFFECTIVE:
                    # destinations
                    # will always seek treatment if symptomatic after treatment failure
                    # will never seek treatment if asymptomatic after treatment failure
                    j = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)
                    if s == SympStat.SYMP.value:
                        dest = ifs_re_tx[j]
                    else:
                        dest = Is[j]

                    ifs_tx_outcome[i] = ChanceNode(
                        name=name,
                        destination_compartments=[dest, dest],
                        probability_params=Constant(1))

    # if susceptible to both CIP and TET
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'If Tx-CIP if test susceptible to CIP and TET | ' \
                   + covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

            dest_cip = ifs_tx_outcome[covert_symp_resist_antibio.get_row_index(
                symp_state=s, rest_profile=p, antibiotic=AB.CIP.value)]
            dest_tet = ifs_tx_outcome[covert_symp_resist_antibio.get_row_index(
                symp_state=s, rest_profile=p, antibiotic=AB.TET.value)]

            ifs_CIP_or_TET[i] = ChanceNode(
                        name=name,
                        destination_compartments=[dest_cip, dest_tet],
                        probability_params=params.probTxCIPIfSuspToCIPAndTET)

    # result of rapid susceptibility test for TET
    # after positive result for rapid CIP test (susceptibility to CIP)
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'Result of TET test after CIP test returned positive | ' \
                   + covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

            # if TET test result is positive (susceptibility to TET)
            # hence the patient can be treated with either CIP or TET
            dest_pos_result = ifs_CIP_or_TET[i]
            # if TET test results is negative (reduced susceptibility to TET)
            # hence the patient will be treated with CIP
            dest_neg_result = ifs_tx_outcome[covert_symp_resist_antibio.get_row_index(
                symp_state=s, rest_profile=p, antibiotic=AB.CIP.value)]

            ifs_rapid_TET_outcome_after_susp_CIP[i] = ChanceNode(
                        name=name,
                        destination_compartments=[dest_pos_result, dest_neg_result],
                        probability_params=params.posTETTest[p])

    # result of rapid susceptibility test for TET
    # after negative result for rapid CIP test (reduced susceptibility to CIP)
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'Result of TET test after CIP test returned negative | ' \
                   + covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

            # if TET test result is positive (susceptibility to TET)
            # hence the patient will be treated with TET
            dest_pos_result = ifs_tx_outcome[covert_symp_resist_antibio.get_row_index(
                symp_state=s, rest_profile=p, antibiotic=AB.TET.value)]
            # if TET test results is negative (reduced susceptibility to TET)
            # hence the patient can be treated with either CFX
            dest_neg_result = ifs_tx_outcome[covert_symp_resist_antibio.get_row_index(
                symp_state=s, rest_profile=p, antibiotic=AB.CFX.value)]

            ifs_rapid_TET_outcome_after_reduced_susp_CIP[i] = ChanceNode(
                name=name,
                destination_compartments=[dest_pos_result, dest_neg_result],
                probability_params=params.posTETTest[p])

    # result of rapid susceptibility test for CIP
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'Result of CIP test | ' + covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

            # if CIP test result is positive (susceptibility to CIP)
            dest_pos_result = ifs_rapid_TET_outcome_after_susp_CIP[i]
            # if CIP test results is negative (reduced susceptibility to CIP)
            dest_neg_result = ifs_rapid_TET_outcome_after_reduced_susp_CIP[i]

            ifs_rapid_CIP_outcome[i] = ChanceNode(
                name=name,
                destination_compartments=[dest_pos_result, dest_neg_result],
                probability_params=params.posCIPTest[p])

    # if will receive a rapid test
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):
            name = 'If will receive a rapid test | ' \
                   + covert_symp_resist.get_str_symp_susp(symp_state=s, rest_profile=p)
            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)

            # if will receive a rapid test
            dest_yes = ifs_rapid_CIP_outcome[i]
            # if will not receive a rapid test, then will be treated with CFX
            dest_no = ifs_tx_outcome[covert_symp_resist_antibio.get_row_index(
                symp_state=s, rest_profile=p, antibiotic=AB.CFX.value)]

            ifs_will_receive_rapid_test[i] = ChanceNode(
                name=name,
                destination_compartments=[dest_yes, dest_no],
                probability_params=params.probRapidTest)

            # if symptomatic
            if s == SympStat.SYMP.value:
                counting_symp.append(ifs_will_receive_rapid_test[i])
            # by resistance profile
            counting_rest_to[p].append(ifs_will_receive_rapid_test[i])

    # ------------- compartment histories ---------------
    # set up prevalence, incidence, and cumulative incidence to collect
    if sets.ifCollectTrajsOfCompartments:
        S.setup_history(collect_prev=True)
        for i in Is:
            i.setup_history(collect_prev=True)
        for f in Fs:
            f.setup_history(collect_prev=True)
        for r in ifs_will_receive_rapid_test:
            r.setup_history(collect_incd=True)
        for r in ifs_rapid_CIP_outcome:
            r.setup_history(collect_incd=True)
        for r in ifs_rapid_TET_outcome_after_susp_CIP:
            r.setup_history(collect_incd=True)
        for r in ifs_rapid_TET_outcome_after_reduced_susp_CIP:
            r.setup_history(collect_incd=True)
        for r in ifs_CIP_or_TET:
            r.setup_history(collect_incd=True)
        for r in ifs_tx_outcome:
            r.setup_history(collect_incd=True)
        for r in ifs_symp_from_emerg_rest:
            r.setup_history(collect_incd=True)
        for r in ifs_re_tx:
            r.setup_history(collect_incd=True)
        for r in ifs_resist_after_re_tx_cfx:
            if r is not None:
                r.setup_history(collect_incd=True)

    for a in range(len(AB)):
        counting_tx_success_by_ab[a].setup_history(collect_incd=True)
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
                           compartments=ifs_will_receive_rapid_test)
    gono_rate = RatioTimeSeries(name='Rate of gonorrhea cases',
                                numerator_sum_time_series=n_cases,
                                denominator_sum_time_series=pop_size,
                                if_surveyed=True,
                                collect_stat_after_warm_up=True)

    # % cases symptomatic
    n_cases_sympt = SumIncidence(name='New cases symptomatic',
                                 compartments=counting_symp)
    perc_cases_sympt = RatioTimeSeries(name='Proportion of cases symptomatic',
                                       numerator_sum_time_series=n_cases_sympt,
                                       denominator_sum_time_series=n_cases,
                                       if_surveyed=True)

    # cases by resistance profile
    n_cases_by_resistance_profile = []
    perc_cases_by_resistance_profile = []
    for p in range(n_rest_profiles):
        n_resistant_cases = SumIncidence(name='Cases resistant to ' + REST_PROFILES[p],
                                         compartments=counting_rest_to[p])
        perc_cases_resistant = RatioTimeSeries(
            name='Proportion of cases resistant to ' + REST_PROFILES[p],
            numerator_sum_time_series=n_resistant_cases,
            denominator_sum_time_series=n_cases,
            if_surveyed=True,
            survey_size_param=params.surveySize)

        n_cases_by_resistance_profile.append(n_resistant_cases)
        perc_cases_by_resistance_profile.append(perc_cases_resistant)

    # cases by resistance to CFX
    n_cases_CFX_R = SumIncidence(
        name='Cases CFX-R, CIP_CFX-R, TET_CFX-R, or CIP_TET_CFX-R',
        compartments=counting_rest_to[RestProfile.CFX.value] +
                     counting_rest_to[RestProfile.CIP_CFX.value] +
                     counting_rest_to[RestProfile.TET.value] +
                     counting_rest_to[RestProfile.CIP_TET_CFX.value])
    perc_cases_CFX_R = RatioTimeSeries(
            name='Proportion of cases CFX-R or CFX+PEN-R',
            numerator_sum_time_series=n_cases_CFX_R,
            denominator_sum_time_series=n_cases,
            if_surveyed=True,
            survey_size_param=params.surveySize)

    # treated with any antibiotics
    n_treated = SumIncidence(
        name='Cases treated',
        compartments=[counting_success_CIP_TET_CFX, counting_tx_M])
    n_treated_CIP_PEN_CFX = SumIncidence(name='Treated with CIP, TET, or CFX',
                                         compartments=[counting_success_CIP_TET_CFX])
    perc_treated_with_CIP_PEN_CFX = RatioTimeSeries(
        name='Proportion of cases treated with CIP, TET, or CFX',
        numerator_sum_time_series=n_treated_CIP_PEN_CFX,
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
    first_line_tx = InterventionAffectingEvents(name='1st line therapy')

    # ------------- attach epidemic events ---------------
    # attached epidemic events to compartments
    # add events to S
    for p in range(n_rest_profiles):
        S.add_event(EpiDepEvent(
                name='Infection | ' + REST_PROFILES[p],
                destination=ifs_symp_from_S[p],
                generating_pathogen=p))

    # add events to infection compartments
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):

            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)
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
                destination=ifs_will_receive_rapid_test[i],
                interv_to_activate=first_line_tx))
            # Is[i].add_event(EpiIndepEvent(
            #     name='Screening then M| ' + compart_name,
            #     rate_param=params.rateScreened,
            #     destination=ifs_counting_tx_M[i],
            #     interv_to_activate=first_line_tx_with_M))
            # seeking treatment
            if s == SympStat.SYMP.value:
                Is[i].add_event(EpiIndepEvent(
                    name='Seeking treatment | ' + compart_name,
                    rate_param=params.rateTreatment,
                    destination=ifs_will_receive_rapid_test[i],
                    interv_to_activate=first_line_tx))
                # Is[i].add_event(EpiIndepEvent(
                #     name='Seeking treatment then M | ' + compart_name,
                #     rate_param=params.rateTreatment,
                #     destination=ifs_counting_tx_M[i],
                #     interv_to_activate=first_line_tx_with_M))

    # add events to infection compartments after treatment failure
    for s in range(n_symp_states):
        for p in range(n_rest_profiles):

            i = covert_symp_resist.get_row_index(symp_state=s, rest_profile=p)
            compart_name = Fs[i].name

            # receive CFX if susceptible to CFX and otherwise, will receive M
            # treatment outcome from receiving CFX
            next_p, reason_for_failure = get_profile_after_resit_or_failure(
                rest_profile=p, antibiotic=AB.CFX)

            # if susceptible to CFX, then will check if resistance could be developed
            if reason_for_failure == TreatmentOutcome.RESISTANCE:
                dest = ifs_resist_after_re_tx_cfx[p]
            # if resistant to CFX, then will receive M
            elif reason_for_failure == TreatmentOutcome.INEFFECTIVE:
                dest = counting_tx_M

            # treatment
            Fs[i].add_event(EpiIndepEvent(
                name='Re-Tx | ' + compart_name,
                rate_param=params.rateRetreatment,
                destination=dest
            ))

    # ------------- populate the model ---------------
    # populate the model
    chance_nodes = ifs_symp_from_S \
                   + ifs_will_receive_rapid_test \
                   + ifs_rapid_CIP_outcome \
                   + ifs_rapid_TET_outcome_after_susp_CIP \
                   + ifs_rapid_TET_outcome_after_reduced_susp_CIP \
                   + ifs_CIP_or_TET + ifs_tx_outcome \
                   + ifs_re_tx + ifs_symp_from_emerg_rest \
                   + counting_tx_success_by_ab\
                   + ifs_resist_after_re_tx_cfx \
                   + [counting_success_CIP_TET_CFX, counting_tx_M]

    list_of_sum_time_series = [pop_size, n_infected, n_cases, n_cases_CFX_R,
                               n_cases_sympt, n_treated, n_treated_CIP_PEN_CFX]
    list_of_sum_time_series.extend(n_cases_by_resistance_profile)

    list_of_ratio_time_series = [prevalence, gono_rate, perc_cases_CFX_R,
                                 perc_cases_sympt, perc_treated_with_CIP_PEN_CFX]
    list_of_ratio_time_series.extend(perc_cases_by_resistance_profile)

    model.populate(compartments=all_comparts,
                   chance_nodes=chance_nodes,
                   list_of_sum_time_series=list_of_sum_time_series,
                   list_of_ratio_time_series=list_of_ratio_time_series,
                   interventions=[first_line_tx],
                   parameters=params)
