from apace.Inputs import ModelSettings
from definitions import get_survey_size, SIM_DURATION, END_OF_WARM_UP, END_OF_CALIB
from model.Data import Prevalence, GonorrheaRate, PercSymptomatic


class GonoSettings(ModelSettings):
    """ settings of the gonorrhea model """

    def __init__(self, if_calibrating=False, collect_traj_of_comparts=True):
        """
        :param if_calibrating: (bool) if calibrating the model
        :param collect_traj_of_comparts: (bool) if collect the trajectories of all compartments
        """

        ModelSettings.__init__(self)

        # model settings
        self.deltaT = 7 / 364  # 1 day
        self.simulationDuration = SIM_DURATION  # years
        self.endOfWarmUpPeriod = END_OF_WARM_UP
        self.simulationOutputPeriod = 1  # simulation output period
        self.observationPeriod = 1
        self.storeParameterValues = True
        self.ifCollectTrajsOfCompartments = collect_traj_of_comparts
        self.exportCalibrationTrajs = False  # if export calibration trajectories

        # calibration settings
        self.calcLikelihood = if_calibrating

        # projection period
        self.storeProjectedOutcomes = True

        # sensitivity and specificity
        self.sensCIP = 0
        self.specCIP = 1
        self.sensTET = 0
        self.specTET = 1

        # probability of receiving a rapid test
        self.probRapidTest = 0 if if_calibrating else 1

        # probability of receiving CIP if someone is susceptible to both CIP and TET
        self.probTxCIPIfSuspToCIPAndTET = 1

        # calibration targets
        if self.calcLikelihood:
            self.simulationDuration = END_OF_CALIB
            self.prevMean = []
            self.prevN = []
            self.gonoRateMean = []
            self.gonoRateN = []
            self.percSympMean = []
            self.percSympN = []

            n_periods_to_keep_constant = 5
            for i in range(len(Prevalence)+n_periods_to_keep_constant-1):

                j = min(i, 1)

                # mean and N of prevalence estimate
                if Prevalence[j][1] is None:
                    self.prevMean.append(None)
                    self.prevN.append(None)
                else:
                    self.prevMean.append(Prevalence[j][1] * 0.01)
                    self.prevN.append(get_survey_size(mean=Prevalence[j][1],
                                                      l=Prevalence[j][2],
                                                      u=Prevalence[j][3],
                                                      multiplier=0.01))
                # mean and N of rate estimate
                if GonorrheaRate[j][1] is None:
                    self.gonoRateMean.append(None)
                    self.gonoRateN.append(None)
                else:
                    self.gonoRateMean.append(GonorrheaRate[j][1] * 0.00001)
                    self.gonoRateN.append(get_survey_size(mean=GonorrheaRate[j][1],
                                                          l=GonorrheaRate[j][1]*0.9,
                                                          u=GonorrheaRate[j][1]*1.1,
                                                          multiplier=0.00001))
                # mean and N of the estimate for the percentage of cases symptomatic
                if PercSymptomatic[j][1] is None:
                    self.percSympMean.append(None)
                    self.percSympN.append(None)
                else:
                    self.percSympMean.append(PercSymptomatic[j][1] * 0.01)
                    self.percSympN.append(get_survey_size(mean=PercSymptomatic[j][1],
                                                          l=PercSymptomatic[j][2],
                                                          u=PercSymptomatic[j][3],
                                                          multiplier=0.01))

    def update_settings(self, sens, spec, prob_rapid_test=0):
        """
        updates certain model parameters and settings
        :param sens: (float) sensitivity of the rapid test for CIP and TET susceptibility
        :param spec: (float) specificity of the rapid test for CIP and TET susceptibility
        :param prob_rapid_test: (float) probability of receiving a rapid test
        """

        self.sensCIP = 0
        self.specCIP = 0
        self.sensTET = 0
        self.specTET = 1
        self.probRapidTest = prob_rapid_test
