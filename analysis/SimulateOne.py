from apace.Epidemic import EpiModel
from model.Model import build_model
from model.ModelSettings import get_model_settings
from model.PlotTrajs import plot

# get model settings
sets = get_model_settings()
# make an (empty) epidemic model
model = EpiModel(id=0, settings=sets)
# populate the SIR model
build_model(model)

# simulate
model.simulate()
# print trajectories
model.export_trajectories(delete_existing_files=True)

# plot trajectories
plot(prev_multiplier=1,    # to show weeks on the x-axis of prevalence data
     incd_multiplier=1*sets.simulationOutputPeriod,     # to show weeks on the x-axis of incidence data
     filename='onetraj'
     )
