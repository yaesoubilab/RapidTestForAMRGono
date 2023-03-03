from apacepy.inputs import ModelSettings
from deampy.in_out_functions import read_csv_rows

from definitions import get_survey_size, SIM_DURATION, END_OF_WARM_UP, END_OF_CALIB, get_scenario_name
from model.data import Prevalence, GonorrheaRate, PercSymptomatic


class GonoSettings(ModelSettings):
    """ settings of the gonorrhea model """

    def __init__(self, if_calibrating=False, collect_traj_of_comparts=True,
                 if_m_available_for_1st_tx=False, sim_duration=None, calibration_seed=None, if_wider_prior=False):
        """
        :param if_calibrating: (bool) if calibrating the model
        :param collect_traj_of_comparts: (bool) if collect the trajectories of all compartments
        :param if_m_available_for_1st_tx: (bool) if M is available for 1st_tx
        :param sim_duration: (float) simulation duration (for sensitivity analysis)
        :param calibration_seed: (int) calibration seed (for sensitivity analysis)
        :param if_wider_prior: (bool) set to True for using wider prior distribution (for sensitivity analysis)
        """

        ModelSettings.__init__(self)

        # model settings
        self.deltaT = 7 / 364  # 1 day
        if sim_duration is None:
            self.simulationDuration = SIM_DURATION  # years
        else:
            self.simulationDuration = sim_duration

        self.endOfWarmUpPeriod = END_OF_WARM_UP
        self.simulationOutputPeriod = 1  # simulation output period
        self.observationPeriod = 1
        self.storeParameterValues = True
        self.ifCollectTrajsOfCompartments = collect_traj_of_comparts
        self.exportCalibrationTrajs = False  # if export calibration trajectories

        # calibration settings
        self.calcLikelihood = if_calibrating
        self.calibSeed = calibration_seed

        # projection period
        self.storeProjectedOutcomes = True

        # sensitivity and specificity
        self.sensCIP = None
        self.specCIP = None
        self.sensTET = None
        self.specTET = None

        # probability of receiving a rapid test
        self.probRapidTest = 0 if if_calibrating else 1

        # adjusting the transmission parameter
        self.transmissionFactor = 1

        # if M is available for the 1st-Tx
        self.ifMAvailableFor1stTx = if_m_available_for_1st_tx
        self.switchThreshold = 0.05

        # folders
        scenario_name = get_scenario_name(if_m_available=self.ifMAvailableFor1stTx,
                                          sim_duration=self.simulationDuration,
                                          calibration_seed=self.calibSeed,
                                          if_wider_priors=if_wider_prior)
        calib_scenario = get_scenario_name(if_m_available=True,
                                           sim_duration=None,
                                           calibration_seed=self.calibSeed,
                                           if_wider_priors=if_wider_prior)

        self.folderToSaveTrajs = 'outputs/{}/trajectories'.format(scenario_name)
        self.folderToSaveSummary = 'outputs/{}/summary'.format(scenario_name)
        self.folderToSaveScenarioAnalysis = 'outputs/{}/scenarios'.format(scenario_name)
        self.folderToSaveCalibrationResults = 'outputs/{}/calibration'.format(calib_scenario)

        # probability of receiving CIP if someone is susceptible to both CIP and TET
        self.probTxCIPIfSuspToCIPAndTET = 0.5

        # calibration targets
        if if_calibrating:
            self.simulationDuration = END_OF_CALIB
            self.ifWiderPrior = if_wider_prior
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

    def update_settings(self, cip_sens, cip_spec, tet_sens, tet_spec, prob_rapid_test=0, transmission_factor=1):
        """
        updates certain model parameters and settings
        :param cip_sens: (float) sensitivity of the rapid test for CIP susceptibility
        :param cip_spec: (float) specificity of the rapid test for CIP susceptibility
        :param tet_sens: (float) sensitivity of the rapid test for TET susceptibility
        :param tet_spec: (float) specificity of the rapid test for TET susceptibility
        :param prob_rapid_test: (float) probability of receiving a rapid test
        :param transmission_factor: (float) to increase or decrease the transmission parameter for sensitivity analysis
        """

        self.sensCIP = cip_sens
        self.specCIP = cip_spec
        self.sensTET = tet_sens
        self.specTET = tet_spec
        self.probRapidTest = prob_rapid_test
        self.transmissionFactor = transmission_factor

    @staticmethod
    def get_list_mean_ci_of_resistance(file_name):

        data = read_csv_rows(file_name=file_name, if_ignore_first_row=False, if_convert_float=True)

        table = []
        for row in data:
            table.append([row[1], [row[2], row[3]]])

        return table
