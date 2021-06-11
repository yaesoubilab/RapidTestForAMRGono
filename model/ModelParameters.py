from SimPy.Parameters import Constant, Equal, Inverse, Product, Division, \
    OneMinus, Uniform, TenToPower, LinearCombination
from apace.Inputs import EpiParameters


class Parameters(EpiParameters):
    """ class to contain the parameters of an SIS model with resistance """
    def __init__(self, model_sets):
        """
        :param model_sets: (ModelSettings)
        """

        EpiParameters.__init__(self)

        one_over_364 = 1/364

        self.popSize = Constant(1000000)
        self.prevI0 = Uniform(0.03, 0.06)
        self.percI0Sym = Uniform(0.0, 0.25)

        self.percI0ResPEN = Uniform(0.0, 0.04)
        self.percI0ResCFX_PEN = Uniform(0.0, 0.005)

        self.transm = Uniform(2, 4)
        self.ratioInfPEN = Uniform(0.8, 1)
        self.ratioInfPEN_CFX = Uniform(0.8, 1)
        self.exponProbResPEN = Uniform(-5, -3)
        self.exponProbResCFX = Uniform(-5, -3)

        self.probSym = Uniform(0.2, 0.8)
        self.tToNaturalRecovery = Uniform(1/12, 5)
        self.tToScreened = Uniform(0.5, 5)
        self.tToTreatment = Uniform(1*one_over_364, 14*one_over_364)
        self.tToRetreatment = Uniform(1 * one_over_364, 14 * one_over_364)
        self.annulSurveySize = Constant(value=1000)

        # calculate dependent parameters
        self.prevS0 = None
        self.prevI0Asym = None
        self.prevI0Res = None
        self.prevI0Sus = None
        self.probResPEN = None
        self.probResCFX = None
        self.rateNaturalRecovery = None
        self.rateScreened = None
        self.rateTreatment = None
        self.rateRetreatment = None
        self.surveySize = None

        self.calculate_dependent_params(model_sets)
        # build the dictionary of parameters
        self.build_dict_of_params()

    def calculate_dependent_params(self, model_sets):

        self.prevS0 = OneMinus(par=self.prevI0)
        self.prevI0Asym = OneMinus(par=self.percI0Sym)
        self.prevI0Res = LinearCombination(parameters=[self.percI0ResPEN, self.percI0ResCFX_PEN])
        self.prevI0Sus = OneMinus(par=self.prevI0Res)

        self.probResPEN = TenToPower(self.exponProbResPEN)
        self.probResCFX = TenToPower(self.exponProbResCFX)

        self.rateNaturalRecovery = Inverse(self.tToNaturalRecovery)
        self.rateScreened = Inverse(self.tToScreened)
        self.rateTreatment = Inverse(self.tToTreatment)
        self.rateRetreatment = Inverse(self.tToRetreatment)

        self.surveySize = Division(self.annulSurveySize, model_sets.observationPeriod)


        self.dictOfParams['Infectivity of susceptible strain'] = Equal(
            par=self.dictOfParams['Transmission parameter'])
        self.dictOfParams['Infectivity of resistant strain'] = Product(
            parameters=[self.dictOfParams['Infectivity of susceptible strain'],
                        self.dictOfParams['Relative infectivity of resistant strain']
                        ])

    def build_dict_of_params(self):
        self.dictOfParams = dict(
            {'Pop size': self.popSize,
             'Initial prevalence': self.prevI0,
             'Initial % I symptomatic': self.percI0Sym,
             'Initial % I PEN-R': self.percI0ResPEN,
             'Initial % I PEN-CFX-R': self.percI0ResCFX_PEN,
             # ----
             'Transmission parameter': self.transm,
             'Relative infectivity of PEN-R strain': self.ratioInfPEN,
             'Relative infectivity of PEN-CFX-R strain': self.ratioInfPEN_CFX,
             # ----
             'Exponent of prob of resistance to PEN': self.exponProbResPEN,
             'Exponent of prob of resistance to CFX': self.exponProbResCFX,
             # ----
             'Prob symptomatic': self.probSym,
             'Time until natural recovery': self.tToNaturalRecovery,
             'Time until screened': self.tToScreened,
             'Time until seeking treatment (symptomatic)': self.tToTreatment,
             'Time until seeking retreatment (symptomatic)': self.tToRetreatment,
             'Annual survey size': self.annulSurveySize
             })

        self.dictOfParams['Initial % Susceptible'] = self.prevS0
        self.dictOfParams['Initial % I asymptomatic'] = self.prevI0Asym
        self.dictOfParams['Initial % I resistant'] = self.prevI0Res
        self.dictOfParams['Initial % I susceptible'] = self.prevI0Sus

        self.dictOfParams['Prob of emergence of PEN-R'] = self.probResPEN
        self.dictOfParams['Prob of emergence of CFX-R'] = self.probResCFX

        self.dictOfParams['Rate of natural recovery'] = self.rateNaturalRecovery
        self.dictOfParams['Rate of screening'] = self.rateScreened
        self.dictOfParams['Rate of seeking treatment'] = self.rateTreatment
        self.dictOfParams['Rate of seeking retreatment'] = self.rateRetreatment

        self.dictOfParams['Survey size (over observation periods)'] = Division(
            par_numerator=self.dictOfParams['Annual survey size'],
            par_denominator=Constant(value=model_sets.observationPeriod))
        self.dictOfParams['Infectivity of susceptible strain'] = Equal(
            par=self.dictOfParams['Transmission parameter'])
        self.dictOfParams['Infectivity of resistant strain'] = Product(
            parameters=[self.dictOfParams['Infectivity of susceptible strain'],
                        self.dictOfParams['Relative infectivity of resistant strain']
                        ])

        # ---------------

        self.dictOfParams['Size of S'] = Product(
            parameters=[self.dictOfParams['Pop size'],
                        self.dictOfParams['Initial % Susceptible']])
        self.dictOfParams['Size of I'] = Product(
            parameters=[self.dictOfParams['Pop size'],
                        self.dictOfParams['Initial prevalence']])

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





