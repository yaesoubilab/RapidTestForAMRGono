from apacepy.inputs import EpiParameters
from deampy.parameters import Constant, Inverse, Product, OneMinus, Uniform, Equal, \
    TenToPower, TimeDependentStepWise, Dirichlet, AnOutcomeOfAMultiVariateDist, \
    ValuesOfParams, TimeDependentSigmoid, Beta

from definitions import RestProfile, AB, SympStat, REST_PROFILES, END_OF_WARM_UP, ConvertSympAndResitAndAntiBio, \
    CIP_SENS_DIST, CIP_SPEC_DIST, TET_SENS_DIST, TET_SPEC_DIST


class Parameters(EpiParameters):
    """ class to contain the parameters of an SIS model with resistance """
    def __init__(self, model_sets):
        """
        :param model_sets: (ModelSettings)
        """

        EpiParameters.__init__(self)

        one_over_364 = 1/364
        self.precIBySymp = [None] * len(SympStat)
        self.percIByRestProfile = [None] * (len(RestProfile) - 1)
        self.ratioInf = [None] * len(RestProfile)
        self.fitnessFMins = [None] * len(RestProfile)
        self.fitnessBs = [None] * len(RestProfile)
        self.fitnessTMids = [None] * len(RestProfile)
        self.exponProbRes = [None] * len(AB)

        # rapid test characteristics
        # beta distributions for sensitivity and specificity
        self.sensCIPBeta = Beta(
            mean=CIP_SENS_DIST[0], st_dev=CIP_SENS_DIST[1], minimum=CIP_SENS_DIST[2], maximum=CIP_SENS_DIST[3])
        self.specCIPBeta = Beta(
            mean=CIP_SPEC_DIST[0], st_dev=CIP_SPEC_DIST[1], minimum=CIP_SPEC_DIST[2], maximum=CIP_SPEC_DIST[3])
        self.sensTETBeta = Beta(
            mean=TET_SENS_DIST[0], st_dev=TET_SENS_DIST[1], minimum=TET_SENS_DIST[2], maximum=TET_SENS_DIST[3])
        self.specTETBeta = Beta(
            mean=TET_SPEC_DIST[0], st_dev=TET_SPEC_DIST[1], minimum=TET_SPEC_DIST[2], maximum=TET_SPEC_DIST[3])

        # use constant if sensitivity/specificity is not provided otherwise, use beta
        self.sensCIP = Equal(self.sensCIPBeta) if model_sets.sensCIP is None else Constant(model_sets.sensCIP)
        self.specCIP = Equal(self.specCIPBeta) if model_sets.specCIP is None else Constant(model_sets.specCIP)
        self.sensTET = Equal(self.sensTETBeta) if model_sets.sensTET is None else Constant(model_sets.sensTET)
        self.specTET = Equal(self.specTETBeta) if model_sets.specTET is None else Constant(model_sets.specTET)

        # if will receive a rapid test
        self.probRapidTest = TimeDependentStepWise(ts=[END_OF_WARM_UP],
                                                   vs=[model_sets.probRapidTest])
        # probability of receiving CIP if someone is susceptible to both CIP and TET
        self.probTxCIPIfSuspToCIPAndTET = Constant(model_sets.probTxCIPIfSuspToCIPAndTET)

        self.popSize = Constant(2.78*10e6)
        self.annulSurveySize = Constant(value=1500)
        self.prevI0 = Uniform(0.03, 0.06)
        self.precIBySymp[SympStat.SYMP.value] = Uniform(0.0, 0.05)

        # the Dirichlet distribution for the percent of I0 by resistance profile
        # (comes from the Excel file under \data folder)
        self.percIByRestProfileDirichlet = Dirichlet(
            par_ns=[55, 2, 270, 0, 170, 0, 0, 0],
            if_ignore_0s=True)

        # this is for debugging purposes ------------
        self.percIConstant = [None] * 8
        self.percIConstant[RestProfile.SUS.value] = Constant(0.8)
        self.percIConstant[RestProfile.CIP.value] = Constant(0.2)
        self.percIConstant[RestProfile.TET.value] = Constant(0.0)
        self.percIConstant[RestProfile.CRO.value] = Constant(0)
        self.percIConstant[RestProfile.CIP_TET.value] = Constant(0)
        self.percIConstant[RestProfile.CIP_CRO.value] = Constant(0)
        self.percIConstant[RestProfile.TET_CRO.value] = Constant(0)
        self.percIConstant[RestProfile.CIP_TET_CRO.value] = Constant(0)
        self.percIByRestProfileConstant = ValuesOfParams(
            parameters=self.percIConstant)
        # --------------------

        # transmission parameters
        self.transmFactor = Constant(model_sets.transmissionFactor) # to adjust transmission for sensitivity analysis
        self.transmBase = Uniform(0.5, 3)  # baseline infectivity
        self.transm = Uniform(0.5, 3)  # baseline infectivity

        # relative infectivity of resistance profiles to susceptible
        if_base_prior = not model_sets.ifWiderPrior
        self.fitnessFMins[RestProfile.SUS.value] = Constant(1)
        self.fitnessFMins[RestProfile.CIP.value] = Uniform(0.9, 1) if if_base_prior else Uniform(0.8, 1)
        self.fitnessFMins[RestProfile.TET.value] = Uniform(0.9, 1) if if_base_prior else Uniform(0.8, 1)
        self.fitnessFMins[RestProfile.CRO.value] = Uniform(0.9, 1) if if_base_prior else Uniform(0.8, 1)
        self.fitnessFMins[RestProfile.CIP_TET.value] = Uniform(0.8, 1) if if_base_prior else Uniform(0.7, 1)
        self.fitnessFMins[RestProfile.CIP_CRO.value] = Uniform(0.8, 1) if if_base_prior else Uniform(0.7, 1)
        self.fitnessFMins[RestProfile.TET_CRO.value] = Uniform(0.8, 1) if if_base_prior else Uniform(0.7, 1)
        self.fitnessFMins[RestProfile.CIP_TET_CRO.value] = Uniform(0.7, 1) if if_base_prior else Uniform(0.6, 1)

        for p in range(len(RestProfile)):
            self.fitnessBs[p] = Uniform(0.1, 0.5) if if_base_prior else Uniform(0.0, 0.6)
            self.fitnessTMids[p] = Uniform(7, 13) if if_base_prior else Uniform(6, 14)

        # exponent of the probability for the emergence of resistance for a drug
        self.exponProbRes[AB.CIP.value] = Uniform(-5, -3) if if_base_prior else Uniform(-6, -2)
        self.exponProbRes[AB.TET.value] = Uniform(-5, -3) if if_base_prior else Uniform(-6, -2)
        self.exponProbRes[AB.CRO.value] = Uniform(-5, -3) if if_base_prior else Uniform(-6, -2)

        self.probSym = Uniform(0.2, 0.8)  # Constant(0.75)
        self.tToNaturalRecovery = Uniform(1/12, 5)  # Constant(4)
        self.tToScreened = Uniform(0.5, 5)  # Constant(4)
        self.tToTreatment = Uniform(1 * one_over_364, 14 * one_over_364)
        self.tToRetreatment = Uniform(1 * one_over_364, 14 * one_over_364)

        # calculate dependent parameters
        # probability of positive CIP or TET tests
        self.posCIPTest = [None] * len(RestProfile)
        self.posTETTest = [None] * len(RestProfile)

        self.prevS = None
        self.percIRes = None
        self.percISus = None
        self.probResEmerge = [None] * len(AB)
        self.rateNaturalRecovery = None
        self.rateScreened = None
        self.rateTreatment = None
        self.rateRetreatment = None
        self.surveySize = None
        self.infectivityByRestProfile = [None] * len(RestProfile)
        self.sizeS = None
        self.sizeI = None
        self.sizeIBySympAndRest = None

        self.calculate_dependent_params(model_sets)
        self.build_dict_of_params()

    def calculate_dependent_params(self, model_sets):

        # transmission
        self.transm = Product(parameters=[self.transmBase, self.transmFactor])

        # percent of I0 by resistance profile
        self.percIByRestProfile = []
        debugging = False
        if not debugging:
            for p in range(len(RestProfile)):
                self.percIByRestProfile.append(AnOutcomeOfAMultiVariateDist(
                    par_multivariate=self.percIByRestProfileDirichlet,
                    outcome_index=p))
        else:
            for p in range(len(RestProfile)):
                self.percIByRestProfile.append(AnOutcomeOfAMultiVariateDist(
                    par_multivariate=self.percIByRestProfileConstant,
                    outcome_index=p))

        # prob of having a positive result for rapid CIP susceptibility test
        # (i.e. the test detects susceptibility to CIP)
        for i, p in enumerate(RestProfile):
            if p in (RestProfile.SUS, RestProfile.TET, RestProfile.CRO, RestProfile.TET_CRO):
                self.posCIPTest[i] = Equal(par=self.sensCIP)
            else:
                self.posCIPTest[i] = OneMinus(par=self.specCIP)

        # prob of having a positive result for rapid TET susceptibility test
        # (i.e., the test detects susceptibility to TET)
        for i, p in enumerate(RestProfile):
            if p in (RestProfile.SUS, RestProfile.CIP, RestProfile.CRO, RestProfile.CIP_CRO):
                self.posTETTest[i] = Equal(par=self.sensTET)
            else:
                self.posTETTest[i] = OneMinus(par=self.specTET)

        self.surveySize = Product([self.annulSurveySize, Constant(model_sets.observationPeriod)])
        self.prevS = OneMinus(par=self.prevI0)
        self.precIBySymp[SympStat.ASYM.value] = OneMinus(par=self.precIBySymp[SympStat.SYMP.value])

        # probability for the emergence of resistance for a drug
        for p in range(len(AB)):
            self.probResEmerge[p] = TenToPower(self.exponProbRes[p])

        self.rateNaturalRecovery = Inverse(self.tToNaturalRecovery)
        self.rateScreened = Inverse(self.tToScreened)
        self.rateTreatment = Inverse(self.tToTreatment)
        self.rateRetreatment = Inverse(self.tToRetreatment)

        # infectivity by susceptibility profile
        for p in range(len(RestProfile)):
            self.ratioInf[p] = TimeDependentSigmoid(
                par_b=self.fitnessBs[p],
                par_t_min=Constant(0),
                par_t_middle=self.fitnessTMids[p],
                par_min=self.fitnessFMins[p],
                par_max=Constant(1))
            self.infectivityByRestProfile[p] = Product([self.transm, self.ratioInf[p]])

        # size of compartments
        self.sizeS = Product([self.popSize, self.prevS])
        self.sizeI = Product([self.popSize, self.prevI0])
        self.sizeIBySympAndRest = [None] * len(SympStat) * len(RestProfile)

        indexer = ConvertSympAndResitAndAntiBio(n_symp_stats=len(SympStat),
                                                n_rest_profiles=len(RestProfile))
        for s in range(len(SympStat)):
            for p in range(len(RestProfile)):
                i = indexer.get_row_index(symp_state=s, rest_profile=p)
                self.sizeIBySympAndRest[i] = Product(
                    [self.sizeI, self.precIBySymp[s], self.percIByRestProfile[p]])

    def build_dict_of_params(self):

        self.dictOfParams = dict(
            {
             'Sensitivity distribution for CIP': self.sensCIPBeta,
             'Specificity distribution for CIP': self.specCIPBeta,
             'Sensitivity distribution for TET': self.sensTETBeta,
             'Specificity distribution for TET': self.specTETBeta,
             'Sensitivity for CIP': self.sensCIP,
             'Specificity for CIP': self.specCIP,
             'Sensitivity for TET': self.sensTET,
             'Specificity for TET': self.specTET,
             'prob CIP test positive': self.posCIPTest,
             'prob TET test positive': self.posTETTest,
             # ---
             'Prob of receiving a rapid test': self.probRapidTest,
             'Prob Tx-CIP if susceptible to CIP and TET': self.probTxCIPIfSuspToCIPAndTET,
             # ----
             'Pop size': self.popSize,
             'Annual survey size': self.annulSurveySize,
             'Initial prevalence': self.prevI0,
             'Initial % I by symptom states': self.precIBySymp,
             'Dirichlet dist. of % I by resistance profile': self.percIByRestProfileDirichlet,
             'Constant dist. of % I': self.percIConstant,
             'Constant dist. of % I by resistance profile': self.percIByRestProfileConstant,

             'Initial % I by resistance profile': self.percIByRestProfile,
             # ----
             'Transmission parameter-base': self.transmBase,
             'Transmission factor': self.transmFactor,
             'Transmission parameter': self.transm,
             'Fitness-f_min': self.fitnessFMins,
             'Fitness-b': self.fitnessBs,
             'Fitness-t_mid': self.fitnessTMids,
             'Relative infectivity by infectivity profile': self.ratioInf,

             # ----
             'Exponent for the prob of resistance by antibiotics': self.exponProbRes,
             # ----
             'Prob symptomatic': self.probSym,
             'Time until natural recovery': self.tToNaturalRecovery,
             'Time until screened': self.tToScreened,
             'Time until seeking treatment (symptomatic)': self.tToTreatment,
             'Time until seeking retreatment (symptomatic)': self.tToRetreatment,
             })

        self.dictOfParams['Survey size (over observation periods)'] = self.surveySize
        self.dictOfParams['Initial % Susceptible'] = self.prevS
        # self.dictOfParams['Initial % I susceptible to all drugs'] = self.percISus

        self.dictOfParams['Prob of resistance by antibiotics'] = self.probResEmerge

        self.dictOfParams['Rate of natural recovery'] = self.rateNaturalRecovery
        self.dictOfParams['Rate of screening'] = self.rateScreened
        self.dictOfParams['Rate of seeking treatment'] = self.rateTreatment
        self.dictOfParams['Rate of seeking retreatment'] = self.rateRetreatment

        for p in range(len(RestProfile)):
            self.dictOfParams['Infectivity of ' + REST_PROFILES[p]] = self.infectivityByRestProfile[p]

        self.dictOfParams['Size of S'] = self.sizeS
        self.dictOfParams['Size of I'] = self.sizeI
        self.dictOfParams['Size of I by Symp/Rest'] = self.sizeIBySympAndRest






