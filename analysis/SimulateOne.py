from apace.Epidemic import EpiModel
from model.Model import build_model
from model.ModelSettings import GonoSettings
from model.Plots import plot

# get model settings
sets = GonoSettings()
sets.update_settings(sens=1, spec=0, prob_rapid_test=1)  # base: (0, 1, 0)

# make an (empty) epidemic model
model = EpiModel(id=1, settings=sets)
# populate the SIR model
build_model(model)

# simulate
model.simulate(seed=1958936649)
# print trajectories
model.export_trajectories(delete_existing_files=True)

# plot trajectories
plot(prev_multiplier=1,    # to show weeks on the x-axis of prevalence data
     incd_multiplier=1*sets.simulationOutputPeriod,     # to show weeks on the x-axis of incidence data
     filename='onetraj'
     )
