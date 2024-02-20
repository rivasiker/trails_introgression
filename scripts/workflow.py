from gwf import Workflow
from os.path import exists
import numpy as np

gwf = Workflow()

mu = 2e-8
g = 25
N_AB = 10_000*2
N_BC = 5_000*2
N_ABC = 10_000*2
N_ref = N_ABC
t_1 = 62_000/g
t_2 = 500_000/g-t_1
t_3 = 6_000_000/g-t_2
t_A = t_1
t_B = t_1
t_C = t_1+t_2
r = 1e-8

t_m = t_1-54_000/g
m = 0.05

dct = {1:8, 3:57, 5:121}

tot_lst = []
t_A = t_1
t_B = t_1
t_C = t_1+t_2
 


model = 'nelder_mead_100Mb'
algorithm = 'Nelder-Mead'
for n_int_AB in [1]:
    for n_int_ABC in [1, 3, 5]:
        if n_int_ABC == 1: mem = 8
        elif n_int_ABC == 3: mem = 20
        elif n_int_ABC == 5: mem = 40
        for seed in range(1, 6):
            gwf.target('simulate_{}_{}_{}_{}'.format(n_int_AB, n_int_ABC, seed, model),
                inputs=['optimize_introgression.py'],
                outputs=['../results/{}_{}_{}_{}_{}.csv'.format(x, n_int_AB, n_int_ABC, seed, model) for x in ['sim']],
                cores=8,
                memory=f'{mem}g',
                walltime= '7-00:00:00',
                account='Primategenomes') << f"""
            python optimize_introgression_100Mb.py {seed} {t_A} {t_B} {t_C} {t_2} {t_3} {N_AB} {N_BC} {N_ABC} {r} {mu} {n_int_AB} {n_int_ABC} {model} {t_m} {m} {algorithm}
            """









# model = 'introgression_error_model'
# algorithm = 'L-BFGS-B'
# for n_int_AB in [1]:
#     for n_int_ABC in [1, 3, 5]:
#         for seed in range(1, 6):
#             # if exists(f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_{model}.csv'):
#             #     pass
#             # print(f"python optimize_introgression.py {seed} {t_A} {t_B} {t_C} {t_2} {t_3} {N_AB} {N_ABC} {r} {mu} {n_int_AB} {n_int_ABC} {model} {t_m} {m}")
#             tot_lst.append('simulate_{}_{}_{}_{}'.format(n_int_AB, n_int_ABC, seed, model))
#             gwf.target('simulate_{}_{}_{}_{}'.format(n_int_AB, n_int_ABC, seed, model),
#                 inputs=['optimize_introgression.py'],
#                 outputs=['../results/{}_{}_{}_{}_{}.csv'.format("sim", n_int_AB, n_int_ABC, seed, model)],
#                 cores=1 if n_int_ABC == 1 else 8,
#                 memory='{}g'.format(n_int_ABC*4),
#                 walltime= 'UNLIMITED',
#                 account='Primategenomes') << f"""
#             python optimize_introgression.py {seed} {t_A} {t_B} {t_C} {t_2} {t_3} {N_AB} {N_ABC} {r} {mu} {n_int_AB} {n_int_ABC} {model} {t_m} {m} {algorithm}
#             """

# model = 'introgression_error_model_NM'
# algorithm = 'Nelder-Mead'
# for n_int_AB in [1]:
#     for n_int_ABC in [1, 3, 5]:
#         for seed in range(1, 6):
#             gwf.target('simulate_{}_{}_{}_{}'.format(n_int_AB, n_int_ABC, seed, model),
#                 inputs=['optimize_introgression.py'],
#                 outputs=['../results/{}_{}_{}_{}_{}.csv'.format(x, n_int_AB, n_int_ABC, seed, model) for x in ['sim']],
#                 cores=1 if n_int_ABC == 1 else 8,
#                 memory='{}g'.format(n_int_ABC*4),
#                 walltime= 'UNLIMITED',
#                 account='Primategenomes') << f"""
#             python optimize_introgression.py {seed} {t_A} {t_B} {t_C} {t_2} {t_3} {N_AB} {N_ABC} {r} {mu} {n_int_AB} {n_int_ABC} {model} {t_m} {m} {algorithm}
#             """


