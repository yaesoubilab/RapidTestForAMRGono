from apace.Epidemic import EpiModel
from model.Model import build_model
from model.ModelSettings import GonoSettings
from model.Plots import plot_trajectories

# get model settings
sets = GonoSettings()
sets.update_settings(sens=1, spec=1, prob_rapid_test=1, # base: (0, 1, 0)
                     if_m_available_for_1st_tx=True)

# make an (empty) epidemic model
model = EpiModel(id=1, settings=sets)
# populate the SIR model
build_model(model)

# simulate
model.simulate(seed=273773726)
# print trajectories
model.export_trajectories(delete_existing_files=True)

# plot trajectories
plot_trajectories(prev_multiplier=1,  # to show weeks on the x-axis of prevalence data
                  incd_multiplier=1*sets.simulationOutputPeriod,  # to show weeks on the x-axis of incidence data
                  filename='onetraj'
                  )
