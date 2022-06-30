from calibrate import N_OF_RESAMPLES_FOR_PARAM_ESTIMATION
from model.support import estimate_parameters

"""
To estimate posterior distribution of model parameters after calibrating the model
and to plot the the posterior distribution and correlation between parameters 
"""

estimate_parameters(n_of_resamples=N_OF_RESAMPLES_FOR_PARAM_ESTIMATION)
