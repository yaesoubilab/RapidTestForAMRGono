from apacepy.analysis.scenarios import ScenarioDataFrame
from deampy.in_out_functions import read_csv_cols_to_dictionary

from definitions import ROOT_DIR, EFFECT_OUTCOME, COST_OUTCOME
from model.scenario_and_sensitivity_analyses import print_corr

param_csvfile = ROOT_DIR + '/analysis/outputs/with M-25yrs/summary/parameter_values.csv'

# selected parameters to include in partial rank correlation analysis
selected_params = [
    'Sensitivity distribution for CIP',
    'Specificity distribution for CIP',
    'Sensitivity distribution for TET',
    'Specificity distribution for TET',
    'Initial prevalence',
    'Initial % I by symptom states-0',
    'Transmission parameter',
    'Exponent for the prob of resistance by antibiotics-0',
    'Exponent for the prob of resistance by antibiotics-1',
    'Exponent for the prob of resistance by antibiotics-2',
    'Prob symptomatic',
    'Time until natural recovery',
    'Time until screened',
    'Time until seeking treatment (symptomatic)',
    'Fitness-f_min-1',
    'Fitness-f_min-2',
    'Fitness-f_min-3',
    'Fitness-f_min-4',
    'Fitness-f_min-5',
    'Fitness-f_min-6',
    'Fitness-f_min-7',
    'Fitness-b-0',
    'Fitness-b-1',
    'Fitness-b-2',
    'Fitness-b-3',
    'Fitness-b-4',
    'Fitness-b-5',
    'Fitness-b-6',
    'Fitness-b-7',
    'Fitness-t_mid-0',
    'Fitness-t_mid-1',
    'Fitness-t_mid-2',
    'Fitness-t_mid-3',
    'Fitness-t_mid-4',
    'Fitness-t_mid-5',
    'Fitness-t_mid-6',
    'Fitness-t_mid-7',
]

# create a dictionary of all parameters along with their values
dict_params = read_csv_cols_to_dictionary(file_name=param_csvfile, if_convert_float=True)

# create a dictionary of selected parameters along with their values
dict_select_params = {k: dict_params[k] for k in selected_params}

# read the simulation scenarios
scenarios = ScenarioDataFrame(
    csv_file_name=ROOT_DIR+'/analysis/outputs/with M-25yrs/scenarios/simulated_scenarios.csv')

# correlation for lifespan
print_corr(scenarios=scenarios, dict_select_params=dict_select_params,
           outcome=EFFECT_OUTCOME, outcome_label='lifespan')

# correlation for gonorrhea rate
print_corr(scenarios=scenarios, dict_select_params=dict_select_params,
           outcome=COST_OUTCOME, outcome_label='gonorrhea rate')
