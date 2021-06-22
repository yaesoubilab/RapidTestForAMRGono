from apace.Inputs import ModelSettings
from definitions import get_survey_size, SIM_DURATION, END_OF_WARM_UP
from model import Data as D


class GonoSettings(ModelSettings):
    """ settings of the gonorrhea model """

    def __init__(self, if_calibrating=False):
        """
        :param if_calibrating: (bool) if calibrating the model
        """

        ModelSettings.__init__(self)

        # model settings
        self.deltaT = 7 / 364  # 1 day
        self.simulationDuration = SIM_DURATION  # years
        self.endOfWarmUpPeriod = END_OF_WARM_UP
        self.simulationOutputPeriod = 1  # simulation output period
        self.observationPeriod = 1
        self.storeParameterValues = True
        self.ifCollectTrajsOfCompartments = True  # if collect the trajectories of all compartments

        # calibration settings
        self.calcLikelihood = if_calibrating
        self.periodBeforeCalibration = 1
        self.calibPeriod = None # WARM_UP + 5   # years

        # projection period
        self.storeProjectedOutcomes = True

        # sensitivity and specificity
        self.sensitivity = 1
        self.specificity = 0

        # calibration targets
        if self.calcLikelihood:
            self.prevMean = []
            self.prevN = []
            self.gonoRateMean = []
            self.gonoRateN = []
            self.percSympMean = []
            self.percSympN = []
            for i in range(5):
                # mean and N of prevalence estimate
                self.prevMean.append(D.Prevalence[0][1] * 0.01)
                self.prevN.append(get_survey_size(mean=D.Prevalence[0][1],
                                                  l=D.Prevalence[0][2],
                                                  u=D.Prevalence[0][3],
                                                  multiplier=0.01))
                # mean and N of rate estimate
                self.gonoRateMean.append(D.GonorrheaRate[0][1] * 0.00001)
                self.gonoRateN.append(get_survey_size(mean=D.GonorrheaRate[0][1],
                                                      l=D.GonorrheaRate[0][1]*0.9,
                                                      u=D.GonorrheaRate[0][1]*1.1,
                                                      multiplier=0.00001))
                # mean and N of the estimate for the percentage of cases symptomatic
                self.percSympMean.append(D.PercSymptomatic[0][1] * 0.01)
                self.percSympN.append(get_survey_size(mean=D.PercSymptomatic[0][1],
                                                      l=D.PercSymptomatic[0][2],
                                                      u=D.PercSymptomatic[0][3],
                                                      multiplier=0.01))

    def update_settings(self, sensitivity, specificity):

        self.sensitivity = sensitivity
        self.specificity = specificity
