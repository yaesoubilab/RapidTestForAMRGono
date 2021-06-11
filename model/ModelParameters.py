from SimPy.Parameters import Constant, Equal, Inverse, Product, Division, \
    OneMinus, Uniform, TenToPower, LinearCombination
from apace.Inputs import EpiParameters
from definitions import SuspProfile, AB, SymStat, SUSP_PROFILES, SympSuspProfiles


class Parameters(EpiParameters):
    """ class to contain the parameters of an SIS model with resistance """
    def __init__(self, model_sets):
        """
        :param model_sets: (ModelSettings)
        """

        EpiParameters.__init__(self)

        one_over_364 = 1/364
        self.percI0BySuspProfile = [None] * len(SuspProfile)
        self.ratioInf = [None] * len(SuspProfile)
        self.exponProbRes = [None] * len(AB)

        self.popSize = Constant(1000000)
        self.annulSurveySize = Constant(value=1000)
        self.prevI0 = Uniform(0.03, 0.06)
        self.percI0Sym = Uniform(0.0, 0.25)

        # percent of I0 by susceptibility profile
        self.percI0BySuspProfile[SuspProfile.SUS.value] = None  # will calculate later
        self.percI0BySuspProfile[SuspProfile.PEN.value] = Uniform(0.0, 0.04)
        self.percI0BySuspProfile[SuspProfile.PEN_CFX.value] = Uniform(0.0, 0.005)

        # infectivity parameters
        self.transm = Uniform(2, 4)  # baseline infectivity
        self.ratioInf[SuspProfile.SUS.value] = Constant(1)
        self.ratioInf[SuspProfile.PEN.value] = Uniform(0.8, 1)
        self.ratioInf[SuspProfile.PEN_CFX.Value] = Uniform(0.8, 1)

        # exponent of the probability for the emergence of resistance for a drug
        for i in range(len(AB)):
            self.exponProbRes[i] = Uniform(-5, -3)

        self.probSym = Uniform(0.2, 0.8)
        self.tToNaturalRecovery = Uniform(1/12, 5)
        self.tToScreened = Uniform(0.5, 5)
        self.tToTreatment = Uniform(1*one_over_364, 14*one_over_364)
        self.tToRetreatment = Uniform(1 * one_over_364, 14 * one_over_364)

        # calculate dependent parameters
        self.prevS0 = None
        self.prevI0Asym = None
        self.prevI0Res = None
        self.probResEmerge = [None] * len(AB)
        self.rateNaturalRecovery = None
        self.rateScreened = None
        self.rateTreatment = None
        self.rateRetreatment = None
        self.surveySize = None
        self.infectivityBySuspProfile = [None] * len(SuspProfile)
        self.sizeS = None
        self.sizeI = None
        self.sizeIBySympAndSusp = None

        self.calculate_dependent_params(model_sets)
        self.build_dict_of_params()

    def calculate_dependent_params(self, model_sets):

        self.surveySize = Division(self.annulSurveySize, model_sets.observationPeriod)
        self.prevS0 = OneMinus(par=self.prevI0)
        self.prevI0Asym = OneMinus(par=self.percI0Sym)

        # find the prevalence of I0 that are susceptible to all antibiotics
        self.prevI0Res = LinearCombination(parameters=self.percI0BySuspProfile[:-1])
        self.percI0BySuspProfile[SuspProfile.SUS.value] = OneMinus(par=self.prevI0Res)

        # probability for the emergence of resistance for a drug
        for i in range(len(AB)):
            self.probResEmerge[i] = TenToPower(self.exponProbRes[i])

        self.rateNaturalRecovery = Inverse(self.tToNaturalRecovery)
        self.rateScreened = Inverse(self.tToScreened)
        self.rateTreatment = Inverse(self.tToTreatment)
        self.rateRetreatment = Inverse(self.tToRetreatment)

        # infectivity by susceptibility profile
        for i in range(len(SuspProfile)):
            self.infectivityBySuspProfile[i] = Product([self.transm, self.ratioInf[i]])

        # size of compartments
        indexer = SympSuspProfiles(n_symp_stats=len(SymStat), n_susp_profiles=len(SuspProfile))
        self.sizeS = Product(self.popSize, self.prevS0)
        self.sizeI = Product(self.popSize, self.prevI0)
        self.sizeIBySympAndSusp = [None] * indexer.length
        for s in range(len(SymStat)):
            for i in range(len(SuspProfile)-1):
                j = indexer.get_row_index(symp_state=s, susp_profile=i)
                self.sizeIBySympAndSusp[j] = Product([self.sizeI, self.percI0Sym])

    def build_dict_of_params(self):

        self.dictOfParams = dict(
            {'Pop size': self.popSize,
             'Annual survey size': self.annulSurveySize,
             'Initial prevalence': self.prevI0,
             'Initial % I symptomatic': self.percI0Sym,
             'Initial % I by susceptibility profile': self.percI0BySuspProfile,
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
        self.dictOfParams['Initial % I asymptomatic'] = self.prevI0Asym
        self.dictOfParams['Initial % I resistant to any drug'] = self.prevI0Res
        self.dictOfParams['Initial % I by susceptibility profile'] = self.percI0BySuspProfile

        self.dictOfParams['Prob of resistance by antibiotics'] = self.probResEmerge

        self.dictOfParams['Rate of natural recovery'] = self.rateNaturalRecovery
        self.dictOfParams['Rate of screening'] = self.rateScreened
        self.dictOfParams['Rate of seeking treatment'] = self.rateTreatment
        self.dictOfParams['Rate of seeking retreatment'] = self.rateRetreatment

        for i in range(len(SuspProfile)):
            self.dictOfParams['Infectivity of '+SUSP_PROFILES[i]] = self.infectivityBySuspProfile[i]

        self.dictOfParams['Size of S'] = self.sizeS
        self.dictOfParams['Size of I'] = self.sizeI


        self.dictOfParams['Size of I-Sus|Sym'] = Product(
            parameters=[self.dictOfParams['Size of I'],
                        self.dictOfParams['Initial % I symptomatic'],
                        self.dictOfParams['Initial % I susceptible']])
        self.dictOfParams['Size of I-Sus|Asym'] = Product(
            parameters=[self.dictOfParams['Size of I'],
                        self.dictOfParams['Initial % I asymptomatic'],
                        self.dictOfParams['Initial % I susceptible']])
        self.dictOfParams['Size of I-Res|Sym'] = Product(
            parameters=[self.dictOfParams['Size of I'],
                        self.dictOfParams['Initial % I symptomatic'],
                        self.dictOfParams['Initial % I resistant']])
        self.dictOfParams['Size of I-Res|Asym'] = Product(
            parameters=[self.dictOfParams['Size of I'],
                        self.dictOfParams['Initial % I asymptomatic'],
                        self.dictOfParams['Initial % I resistant']])





