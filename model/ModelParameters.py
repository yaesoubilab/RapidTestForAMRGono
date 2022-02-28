from SimPy.Parameters import Constant, Inverse, Product, OneMinus, Uniform, Equal, \
    TenToPower, OneMinusSum
from apace.Inputs import EpiParameters
from definitions import RestProfile, AB, SympStat, REST_PROFILES, ConvertSympAndSuspAndAntiBio


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
        self.exponProbRes = [None] * len(AB)

        # rapid test characteristics
        self.sensCIP = Constant(model_sets.sensCIP)
        self.specCIP = Constant(model_sets.specCIP)
        self.sensTET = Constant(model_sets.sensTET)
        self.specTET = Constant(model_sets.specTET)

        # if will receive a rapid test
        self.probRapidTest = Constant(model_sets.probRapidTest)
        # probability of receiving CIP if someone is susceptible to both CIP and TET
        self.probTxCIPIfSuspToCIPAndTET = Constant(model_sets.probTxCIPIfSuspToCIPAndTET)

        self.popSize = Constant(1000000)
        self.annulSurveySize = Constant(value=1000)
        self.prevI0 = Uniform(0.03, 0.06)
        self.precIBySymp[SympStat.SYMP.value] = Uniform(0.0, 0.05)

        # percent of I0 by resistance profile (comes from the Excel file ResistanceData.xlms)
        self.percIByRestProfile[RestProfile.CIP.value] = Uniform(0.002, 0.007)
        self.percIByRestProfile[RestProfile.TET.value] = Uniform(0.448, 0.672)
        self.percIByRestProfile[RestProfile.CFX.value] = Uniform(0, 0)
        self.percIByRestProfile[RestProfile.CIP_TET.value] = Uniform(0.126, 0.189)
        self.percIByRestProfile[RestProfile.CIP_CFX.value] = Uniform(0, 0)
        self.percIByRestProfile[RestProfile.TET_CFX.value] = Uniform(0, 0)
        self.percIByRestProfile[RestProfile.CIP_TET_CFX.value] = Uniform(0, 0)

        # infectivity parameters
        self.transm = Uniform(0.5, 3)  # baseline infectivity
        # relative infectivity of resistance profiles to susceptible
        self.ratioInf[RestProfile.CIP.value] = Uniform(0.9, 1)
        self.ratioInf[RestProfile.TET.value] = Uniform(0.9, 1)
        self.ratioInf[RestProfile.CFX.value] = Uniform(0.9, 1)
        self.ratioInf[RestProfile.CIP_TET.value] = Uniform(0.8, 1)
        self.ratioInf[RestProfile.CIP_CFX.value] = Uniform(0.8, 1)
        self.ratioInf[RestProfile.TET_CFX.value] = Uniform(0.8, 1)
        self.ratioInf[RestProfile.CIP_TET_CFX.value] = Uniform(0.7, 1)
        self.ratioInf[RestProfile.SUS.value] = Constant(1)

        # exponent of the probability for the emergence of resistance for a drug
        self.exponProbRes[AB.CIP.value] = Uniform(-5, -3)
        self.exponProbRes[AB.TET.value] = Uniform(-5, -3)
        self.exponProbRes[AB.CFX.value] = Uniform(-5, -3)

        self.probSym = Uniform(0.2, 0.8)  # Constant(0.75)
        self.tToNaturalRecovery = Uniform(1/12, 5)  # Constant(4)
        self.tToScreened = Uniform(0.5, 5)  # Constant(4)
        self.tToTreatment = Uniform(1 * one_over_364, 14 * one_over_364)
        self.tToRetreatment = Uniform(1 * one_over_364, 14 * one_over_364)

        # calculate dependent parameters
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

        # prob of having a positive result for rapid CIP susceptibility test
        for i, p in enumerate(RestProfile):
            if p in (RestProfile.SUS, RestProfile.TET, RestProfile.CFX, RestProfile.TET_CFX):
                self.posCIPTest[i] = Equal(par=self.sensCIP)
            else:
                self.posCIPTest[i] = OneMinus(par=self.specCIP)

        # prob of having a positive result for rapid TET susceptibility test
        for i, p in enumerate(RestProfile):
            if p in (RestProfile.SUS, RestProfile.CIP, RestProfile.CFX, RestProfile.CIP_CFX):
                self.posTETTest[i] = Equal(par=self.sensTET)
            else:
                self.posTETTest[i] = OneMinus(par=self.specTET)

        self.surveySize = Product([self.annulSurveySize, Constant(model_sets.observationPeriod)])
        self.prevS = OneMinus(par=self.prevI0)
        self.precIBySymp[SympStat.ASYM.value] = OneMinus(par=self.precIBySymp[SympStat.SYMP.value])

        # find the prevalence of I0 that are susceptible to all antibiotics
        self.percISus = OneMinusSum(parameters=self.percIByRestProfile)

        # probability for the emergence of resistance for a drug
        for p in range(len(AB)):
            self.probResEmerge[p] = TenToPower(self.exponProbRes[p])

        self.rateNaturalRecovery = Inverse(self.tToNaturalRecovery)
        self.rateScreened = Inverse(self.tToScreened)
        self.rateTreatment = Inverse(self.tToTreatment)
        self.rateRetreatment = Inverse(self.tToRetreatment)

        # infectivity by susceptibility profile
        for p in range(len(RestProfile)):
            self.infectivityByRestProfile[p] = Product([self.transm, self.ratioInf[p]])

        # size of compartments
        self.sizeS = Product([self.popSize, self.prevS])
        self.sizeI = Product([self.popSize, self.prevI0])
        self.sizeIBySympAndRest = [None] * len(SympStat) * len(RestProfile)

        indexer = ConvertSympAndSuspAndAntiBio(n_symp_stats=len(SympStat),
                                               n_rest_profiles=len(RestProfile))
        for s in range(len(SympStat)):
            for p in range(len(RestProfile)):
                i = indexer.get_row_index(symp_state=s, rest_profile=p)
                if p == RestProfile.SUS.value:
                    self.sizeIBySympAndRest[i] = Product(
                        [self.sizeI, self.precIBySymp[s], self.percISus])
                else:
                    self.sizeIBySympAndRest[i] = Product(
                        [self.sizeI, self.precIBySymp[s], self.percIByRestProfile[p]])

    def build_dict_of_params(self):

        self.dictOfParams = dict(
            {'Sensitivity for CIP': self.sensCIP,
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
             'Initial % I by resistance profile': self.percIByRestProfile,
             # ----
             'Transmission parameter': self.transm,
             'Relative infectivity by susceptibility profile': self.ratioInf,
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
        self.dictOfParams['Initial % I susceptible to all drugs'] = self.percISus

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






