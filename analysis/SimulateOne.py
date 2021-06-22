from apace.Epidemic import EpiModel
from model.Model import build_model
from model.ModelSettings import GonoSettings
from model.Plots import plot

# get model settings
sets = GonoSettings()
# make an (empty) epidemic model
model = EpiModel(id=1, settings=sets)
# populate the SIR model
build_model(model)

# simulate
model.simulate(seed=484264043)
# print trajectories
model.export_trajectories(delete_existing_files=True)

# plot trajectories
plot(prev_multiplier=1,    # to show weeks on the x-axis of prevalence data
     incd_multiplier=1*sets.simulationOutputPeriod,     # to show weeks on the x-axis of incidence data
     filename='onetraj'
     )
