from apacepy.epidemic import EpiModel

from model.model_settings import GonoSettings
from model.model_structure import build_model
from model.plots import plot_trajectories

# if drug M can be used for 1st line therapy
IF_M_AVAILABLE_FOR_FIRST_TX = True
# sensitivity, specificity, and coverage of the rapid test
# the status quo is (0, 1, 0, 1, 0)
CIP_SENS, CIP_SPEC, TET_SENS, TET_SPEC, COVERAGE = 0, 1, 0, 1, 0

# -------------------
# get model settings
sets = GonoSettings(if_m_available_for_1st_tx=IF_M_AVAILABLE_FOR_FIRST_TX)
sets.update_settings(cip_sens=CIP_SENS, cip_spec=CIP_SPEC,
                     tet_sens=TET_SENS, tet_spec=TET_SPEC,
                     prob_rapid_test=COVERAGE)

# make an (empty) epidemic model
model = EpiModel(id=1, settings=sets)
# populate the model with the gonorrhea model
build_model(model)

# simulate
model.simulate(seed=625657653)
# export trajectories
model.export_trajectories()

# file name
if IF_M_AVAILABLE_FOR_FIRST_TX:
    filename = 'onetraj-with-M'
else:
    filename = 'onetraj-no-M'

# plot trajectories
plot_trajectories(prev_multiplier=1,  # to show weeks on the x-axis of prevalence data
                  incd_multiplier=1*sets.simulationOutputPeriod,  # to show weeks on the x-axis of incidence data
                  dir_of_traj_files=sets.folderToSaveTrajs,
                  filename=filename
                  )
