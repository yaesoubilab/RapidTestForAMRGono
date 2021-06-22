from SimPy.Parameters import Constant, Inverse, Product, OneMinus, Uniform, \
    TenToPower, OneMinusSum
from apace.Inputs import EpiParameters
from definitions import RestProfile, AB, SympStat, REST_PROFILES


class Parameters(EpiParameters):
    """ class to contain the parameters of an SIS model with resistance """
    def __init__(self, model_sets):
        """
        :param model_sets: (ModelSettings)
        """

        EpiParameters.__init__(self)

        one_over_364 = 1/364
        self.precI0BySymp = [None] * len(SympStat)
        self.percI0ByRestProfile = [None] * (len(RestProfile) - 1)
        self.ratioInf = [None] * len(RestProfile)
        self.exponProbRes = [None] * len(AB)

        # rapid test characteristics
        self.sens = Constant(model_sets.sensitivity)
        self.spec = Constant(model_sets.specificity)

        self.popSize = Constant(1000000)
        self.annulSurveySize = Constant(value=1000)
        self.prevI0 = Uniform(0.03, 0.06)
        self.precI0BySymp[SympStat.SYMP.value] = Uniform(0.0, 0.1)

        # percent of I0 by resistance profile
        self.percI0ByRestProfile[RestProfile.PEN.value] = Uniform(0.05, 0.15)
        self.percI0ByRestProfile[RestProfile.CFX.value] = Uniform(0.0, 0.001)
        self.percI0ByRestProfile[RestProfile.PEN_CFX.value] = Uniform(0.0, 0.001)

        # infectivity parameters
        self.transm = Uniform(0.5, 3)  # baseline infectivity
        self.ratioInf[RestProfile.PEN.value] = Uniform(0.8, 1)
        self.ratioInf[RestProfile.CFX.value] = Uniform(0.8, 1)
        self.ratioInf[RestProfile.PEN_CFX.value] = Uniform(0.8, 1)
        self.ratioInf[RestProfile.SUS.value] = Constant(1)

        # exponent of the probability for the emergence of resistance for a drug
        self.exponProbRes[AB.PEN.value] = Uniform(-5, -3)
        self.exponProbRes[AB.CFX.value] = Uniform(-5, -3)

        self.probSym = Uniform(0.2, 0.8)  # Constant(0.75)
        self.tToNaturalRecovery = Uniform(1/12, 5)  # Constant(4)
        self.tToScreened = Uniform(0.5, 5)  # Constant(4)
        self.tToTreatment = Uniform(1 * one_over_364, 14 * one_over_364)
        self.tToRetreatment = Uniform(1 * one_over_364, 14 * one_over_364)

        # calculate dependent parameters
        self.oneMinusSpec = None
        self.prevS0 = None
        self.percI0Res = None
        self.percI0Sus = None
        self.probResEmerge = [None] * len(AB)
        self.rateNaturalRecovery = None
        self.rateScreened = None
        self.rateTreatment = None
        self.rateRetreatment = None
        self.surveySize = None
        self.infectivityBySuspProfile = [None] * len(RestProfile)
        self.sizeS = None
        self.sizeI = None
        self.sizeIBySympAndSusp = None

        self.calculate_dependent_params(model_sets)
        self.build_dict_of_params()

    def calculate_dependent_params(self, model_sets):

        self.oneMinusSpec = OneMinus(par=self.spec)
        self.surveySize = Product([self.annulSurveySize, Constant(model_sets.observationPeriod)])
        self.prevS0 = OneMinus(par=self.prevI0)
        self.precI0BySymp[SympStat.ASYM.value] = OneMinus(par=self.precI0BySymp[SympStat.SYMP.value])

        # find the prevalence of I0 that are susceptible to all antibiotics
        self.percI0Sus = OneMinusSum(parameters=self.percI0ByRestProfile)

        # probability for the emergence of resistance for a drug
        for p in range(len(AB)):
            self.probResEmerge[p] = TenToPower(self.exponProbRes[p])

        self.rateNaturalRecovery = Inverse(self.tToNaturalRecovery)
        self.rateScreened = Inverse(self.tToScreened)
        self.rateTreatment = Inverse(self.tToTreatment)
        self.rateRetreatment = Inverse(self.tToRetreatment)

        # infectivity by susceptibility profile
        for p in range(len(RestProfile)):
            self.infectivityBySuspProfile[p] = Product([self.transm, self.ratioInf[p]])

        # size of compartments
        self.sizeS = Product([self.popSize, self.prevS0])
        self.sizeI = Product([self.popSize, self.prevI0])
        self.sizeIBySympAndSusp = [None] * len(SympStat) * len(RestProfile)
        i = 0
        for s in range(len(SympStat)):
            for p in range(len(RestProfile)):
                if p == RestProfile.SUS.value:
                    self.sizeIBySympAndSusp[i] = Product(
                        [self.sizeI, self.precI0BySymp[s], self.percI0Sus])
                else:
                    self.sizeIBySympAndSusp[i] = Product(
                        [self.sizeI, self.precI0BySymp[s], self.percI0ByRestProfile[p]])
                i += 1

    def build_dict_of_params(self):

        self.dictOfParams = dict(
            {'Sensitivity': self.sens,
             'Specificity': self.spec,
             '1-Specificity': self.oneMinusSpec,
             # ----
             'Pop size': self.popSize,
             'Annual survey size': self.annulSurveySize,
             'Initial prevalence': self.prevI0,
             'Initial % I by symptom states': self.precI0BySymp,
             'Initial % I by resistance profile': self.percI0ByRestProfile,
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
        self.dictOfParams['Initial % Susceptible'] = self.prevS0
        self.dictOfParams['Initial % I susceptible to all drugs'] = self.percI0Sus

        self.dictOfParams['Prob of resistance by antibiotics'] = self.probResEmerge

        self.dictOfParams['Rate of natural recovery'] = self.rateNaturalRecovery
        self.dictOfParams['Rate of screening'] = self.rateScreened
        self.dictOfParams['Rate of seeking treatment'] = self.rateTreatment
        self.dictOfParams['Rate of seeking retreatment'] = self.rateRetreatment

        for p in range(len(RestProfile)):
            self.dictOfParams['Infectivity of ' + REST_PROFILES[p]] = self.infectivityBySuspProfile[p]

        self.dictOfParams['Size of S'] = self.sizeS
        self.dictOfParams['Size of I'] = self.sizeI
        self.dictOfParams['Size of I by Symp/Susp'] = self.sizeIBySympAndSusp






