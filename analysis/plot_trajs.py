from model.plots import plot_trajectories


# plot trajectories
plot_trajectories(prev_multiplier=1,  # to show weeks on the x-axis of prevalence data
                  incd_multiplier=1,  # to show weeks on the x-axis of incidence data
                  obs_prev_multiplier=1,
                  obs_incd_multiplier=1,
                  filename='Calibrated.png')
