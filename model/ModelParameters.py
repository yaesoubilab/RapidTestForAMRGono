from SimPy.Parameters import Constant, Equal, Inverse, Product, Division, OneMinus, Uniform, TenToPower, Surge
from apace.Inputs import EpiParameters


class Parameters(EpiParameters):
    """ class to contain the parameters of an SIS model with resistance """
    def __init__(self, model_sets):
        """
        :param model_sets: (ModelSettings)
        """

        EpiParameters.__init__(self)

        one_over_364 = 1/364
        self.dictOfParams = dict(
            {'Pop size': Constant(value=1000000),
             'Initial prevalence': Uniform(0.03, 0.06),  # Constant(value=0.04),
             'Initial % I symptomatic': Uniform(0.0, 0.25),  # Constant(value=0.1),
             'Initial % I resistant': Uniform(0.0, 0.04),  # Constant(value=0.02),
             # ----
             'Transmission parameter': Uniform(2, 4),  # Constant(value=3),
             'Relative infectivity of resistant strain': Uniform(0.8, 1),  # Constant(value=0.8),
             # ----
             'Prob symptomatic': Uniform(0.2, 0.8),  # Constant(value=0.5),
             'Exponent of prob of resistance': Uniform(-5, -3),
             # ----
             'Time until natural recovery': Uniform(1/12, 5),  # Constant(value=2),
             'Time until screened': Uniform(0.5, 5),  # Constant(value=1),
             'Time until seeking treatment (symptomatic)': Uniform(1*one_over_364, 14*one_over_364),  # Constant(value=14 / 364),
             'Time until seeking retreatment (symptomatic)': Uniform(1*one_over_364, 14*one_over_364),  # Constant(value=14 / 364),
             'Annual survey size': Constant(value=1000)
             })

        # period of disruption
        if model_sets.withDisruption:
            t0 = model_sets.periodBeforeCalibration + model_sets.calibPeriod
            t1 = t0 + model_sets.disruptionPeriod
            self.dictOfParams['(time-dep) Infectivity of susceptible strain'] = Surge(
                par_base=self.dictOfParams['Transmission parameter'],
                par_max_percent_change=-0.5, par_t0=t0, par_t1=t1)
            self.dictOfParams['(time-dep) Time until screened'] = Surge(
                par_base=self.dictOfParams['Time until screened'],
                par_max_percent_change=1, par_t0=t0, par_t1=t1)
            self.dictOfParams['(time-dep) Time until seeking treatment (symptomatic)'] = Surge(
                par_base=self.dictOfParams['Time until seeking treatment (symptomatic)'],
                par_max_percent_change=1, par_t0=t0, par_t1=t1)
            self.dictOfParams['(time-dep) Time until seeking retreatment (symptomatic)'] = Surge(
                par_base=self.dictOfParams['Time until seeking retreatment (symptomatic)'],
                par_max_percent_change=1, par_t0=t0, par_t1=t1)
            self.dictOfParams['(time-dep) Annual survey size'] = Surge(
                par_base=self.dictOfParams['Annual survey size'],
                par_max_percent_change=-0.5, par_t0=t0, par_t1=t1)

        self.dictOfParams['Initial % Susceptible'] = OneMinus(par=self.dictOfParams['Initial prevalence'])
        self.dictOfParams['Initial % I asymptomatic'] = OneMinus(par=self.dictOfParams['Initial % I symptomatic'])
        self.dictOfParams['Initial % I susceptible'] = OneMinus(
            par=self.dictOfParams['Initial % I resistant'])

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

        self.dictOfParams['Prob of resistance'] = TenToPower(
            par=self.dictOfParams['Exponent of prob of resistance'])
        self.dictOfParams['Rate of natural recovery'] = Inverse(
            par=self.dictOfParams['Time until natural recovery'])

        if model_sets.withDisruption:
            self.dictOfParams['Rate of screening'] = Inverse(
                par=self.dictOfParams['(time-dep) Time until screened'])
            self.dictOfParams['Rate of seeking treatment'] = Inverse(
                par=self.dictOfParams['(time-dep) Time until seeking treatment (symptomatic)'])
            self.dictOfParams['Rate of seeking retreatment'] = Inverse(
                par=self.dictOfParams['(time-dep) Time until seeking retreatment (symptomatic)'])
            self.dictOfParams['Survey size (over observation periods)'] = Division(
                par_numerator=self.dictOfParams['(time-dep) Annual survey size'],
                par_denominator=Constant(value=model_sets.observationPeriod))
            self.dictOfParams['Infectivity of susceptible strain'] = Equal(
                par=self.dictOfParams['(time-dep) Infectivity of susceptible strain'])
        else:
            self.dictOfParams['Rate of screening'] = Inverse(
                par=self.dictOfParams['Time until screened'])
            self.dictOfParams['Rate of seeking treatment'] = Inverse(
                par=self.dictOfParams['Time until seeking treatment (symptomatic)'])
            self.dictOfParams['Rate of seeking retreatment'] = Inverse(
                par=self.dictOfParams['Time until seeking retreatment (symptomatic)'])
            self.dictOfParams['Survey size (over observation periods)'] = Division(
                par_numerator=self.dictOfParams['Annual survey size'],
                par_denominator=Constant(value=model_sets.observationPeriod))
            self.dictOfParams['Infectivity of susceptible strain'] = Equal(
                par=self.dictOfParams['Transmission parameter'])

        self.dictOfParams['Infectivity of resistant strain'] = Product(
            parameters=[self.dictOfParams['Infectivity of susceptible strain'],
                        self.dictOfParams['Relative infectivity of resistant strain']
                        ])
