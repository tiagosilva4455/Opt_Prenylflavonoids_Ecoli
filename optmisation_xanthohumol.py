import time
from builtins import print
from mewpy.util.constants import EAConstants
from mewpy import solvers
from mewpy.optimization import EA
from mewpy.optimization import set_default_engine
from mewpy.optimization.evaluation import BPCY, WYIELD, ModificationType
from mewpy.problems import GOUProblem
from mewpy.simulation import get_simulator, SimulationMethod
from mewpy.util.constants import EAConstants
from mewpy.util.io import population_to_csv
from reframed.io.sbml import load_cbmodel
from cobra.io import read_sbml_model, validate_sbml_model
import warnings

def main(model_filename,
         biomass_id='BIOMASS_Ec_iML1515_core_75p37M',
         product_id='R08004',
         env_conditions=None):
    if env_conditions is None:
        return

    warnings.filterwarnings("ignore")

    print('Loading model...')
    model = read_sbml_model(model_filename)
    #model = load_cbmodel(model_filename)
    print('Model loaded!')

    print('Setting model objective...')
    model.objective = biomass_id
    print('Model objective set!')

    #EAConstants.LEVELS.append(1)
    # uses 16 parallel threads
    # EAConstants.NUM_CPUS = 6

    s = solvers.get_default_solver()
    print('Default solver to be used: ' + str(s))

    #model.genes.YLR134W.knock_out()
    #model.genes.YDR380W.knock_out()

    for envcond in env_conditions:
        try:
            startTime=time.strftime('%d-%m-%Y__%H_%M_%S')
            print(startTime)
            print(f'Starting with envcond: {envcond}')
            carbon_source = next(iter(envcond))
            print('Setting evaluation functions...')
            evaluator_1 = WYIELD(biomass_id, product_id)
            evaluator_2 = BPCY(biomass_id,
                               product_id,
                               uptake=carbon_source,
                               method=SimulationMethod.pFBA)
            evaluator_3 = ModificationType()
            print('Evaluation functions set!')

            print('Setting the GOUProblem...')
            solution = {'b0754': 4, 'b2600': 4, 'b1702': 4, 'b2935': 4, 'b0908': 4, 'b0388': 4,'b2329': 4, 'b3281': 4, 'b1693': 4, 'b3389': 4}
            # The reaction up and down regulation optimization problem
            problem = GOUProblem(model,
                                 fevaluation=[evaluator_2,
                                              evaluator_1,
                                              evaluator_3
                                              ],
                                 envcond=envcond,
                                 candidate_min_size=2,
                                 candidate_max_size=8,
                                 partial_solution=solution)
            print('GOUProblem set!')

            # validate whether product is achievable
            res = problem.simulate(objective={product_id:1})
            if res.objective_value < 10E-6:
                continue
            print(res)

            print('Flux simulation complete!')

            model.objective = biomass_id

            print('Starting optimization...')

            ea = EA(problem, max_generations=2)

            ea.run(simplify=False)
            print('Optimization complete!')

            endTime = time.strftime('%d-%m-%Y__%H_%M_%S')
            
            print(endTime)

            print('Writing the results...')
            filename = './results/results_'+startTime+'_'+endTime
            for key, value in envcond.items():
                filename += f'{key}_{value}'
            filename += '.csv'
            df=ea.dataframe()
            df.to_csv(filename)
            
            print('End!')

        except BaseException as ex:
            print('Exception occurred!')
            print(ex)


if __name__ == '__main__':
    PRODUCT_ID = 'R08004'

    model_f = 'edited_e_coli_w_knockouts.xml'

    BIOMASS_ID = 'BIOMASS_Ec_iML1515_core_75p37M'
    GLC = 'EX_glc__D_e'
    #XYL = 'r_1718'
    #GLY = 'r_1808'
    O2 = 'EX_o2_e'
    #AMO = 'r_1654'
    #PHO = 'r_2005'

    # model_f = 'iMM904_chond.xml'
    # #model_f = 'iND750_chond.xml'
    # BIOMASS_ID = 'BIOMASS_SC5_notrace'
    # #BIOMASS_ID = 'R_BIOMASS_SC4_bal'
    # GLC = 'EX_glc__D_e'
    # XYL = 'EX_xyl__D_e'
    # GLY = 'EX_glyc_e'
    # PHO = 'EX_pi_e'
    # O2 = 'EX_o2_e'

    env_conds = [  {GLC: (-10.0, 999999.0)},
                  # {GLC: (-10.0, 999999.0), O2: (-1, 100000.0)},
                  # {GLY: (-20.0, 999999.0), GLC: (-10.0, 999999.0), O2: (-1, 100000.0)},
                  # {GLC: (-10.0, 999999.0), GLY: (-20.0, 999999.0), PHO: (-1.0, 999999.0)},
                  # {GLC: (-10.0, 999999.0), PHO: (-1.0, 999999.0)},
                  # {XYL: (-12.0, 999999.0)},
                  # {XYL: (-12.0, 999999.0), O2: (-1, 100000.0)},
                  # {GLY: (-20.0, 999999.0), O2: (-1, 100000.0)},
                  # {GLY: (-20.0, 999999.0)}
                  ]

    main(model_f, biomass_id=BIOMASS_ID, product_id=PRODUCT_ID, env_conditions=env_conds)
