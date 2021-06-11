from apace.Inputs import ModelSettings
from definitions import get_survey_size
from model import Data as D


def get_model_settings(if_calibrating=False, if_disruption=False):
    """
    :param if_calibrating: (bool) if calibrating the model
    :param if_disruption: (bool) if simulate trajectories under disruptions
    :return:
    """

    # model settings
    settings = ModelSettings()
    settings.deltaT = 7 / 364  # 1 day
    settings.simulationDuration = 10  # years
    settings.simulationOutputPeriod = 1  # simulation output period
    settings.observationPeriod = 1
    settings.storeParameterValues = True

    # calibration settings
    settings.calcLikelihood = if_calibrating
    settings.periodBeforeCalibration = 1
    settings.calibPeriod = 5  # years

    # projection period
    settings.disruptionPeriod = 2  # years
    settings.withDisruption = if_disruption

    # calibration targets
    if settings.calcLikelihood:
        settings.prevMean = []
        settings.prevN = []
        settings.gonoRateMean = []
        settings.gonoRateN = []
        settings.percSympMean = []
        settings.percSympN = []
        for i in range(5):
            # mean and N of prevalence estimate
            settings.prevMean.append(D.Prevalence[0][1] * 0.01)
            settings.prevN.append(get_survey_size(mean=D.Prevalence[0][1],
                                                  l=D.Prevalence[0][2],
                                                  u=D.Prevalence[0][3],
                                                  multiplier=0.01))
            # mean and N of rate estimate
            settings.gonoRateMean.append(D.GonorrheaRate[0][1] * 0.00001)
            settings.gonoRateN.append(get_survey_size(mean=D.GonorrheaRate[0][1],
                                                      l=D.GonorrheaRate[0][1]*0.8,
                                                      u=D.GonorrheaRate[0][1]*1.2,
                                                      multiplier=0.00001))
            # mean and N of the estimate for the percentage of cases symptomatic
            settings.percSympMean.append(D.PercSymptomatic[0][1] * 0.01)
            settings.percSympN.append(get_survey_size(mean=D.PercSymptomatic[0][1],
                                                      l=D.PercSymptomatic[0][2],
                                                      u=D.PercSymptomatic[0][3],
                                                      multiplier=0.01))

    return settings
