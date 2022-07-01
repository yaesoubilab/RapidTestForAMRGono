from model.plots import plot_trajectories


"""
To plot the simulated trajectories available under 
"""

# plot trajectories
plot_trajectories(prev_multiplier=1,  # to show weeks on the x-axis of prevalence data
                  incd_multiplier=1,  # to show weeks on the x-axis of incidence data
                  obs_prev_multiplier=1,
                  obs_incd_multiplier=1,
                  dir_of_trajs='outputs/with-M/trajectories',
                  filename='Calibrated.png')
