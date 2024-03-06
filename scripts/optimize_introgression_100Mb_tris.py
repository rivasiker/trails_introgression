import sys
from trails.optimizer import trans_emiss_calc_introgression
from trails.cutpoints import cutpoints_ABC, cutpoints_AB
import numpy as np
from trails.optimizer import loglik_wrapper_par, write_list, optimizer_introgression
from trails.read_data import get_obs_state_dct
import pandas as pd
import time
import re
import msprime
import sys

####################### Model parameters #######################

print('Model parameters')

n_sites = 10_000_000

seed = int(sys.argv[1])
t_A = float(sys.argv[2])
t_B = float(sys.argv[3])
t_C = float(sys.argv[4])
t_2 = float(sys.argv[5])
t_3 = float(sys.argv[6])
N_AB = float(sys.argv[7])
N_BC = float(sys.argv[8])
N_ABC = float(sys.argv[9])
r = float(sys.argv[10])
mu = float(sys.argv[11])
n_int_AB = int(sys.argv[12])
n_int_ABC = int(sys.argv[13])
model = str(sys.argv[14])

t_m = float(sys.argv[15])
m = float(sys.argv[16])

algorithm = str(sys.argv[17])

t_C_prime = t_C-t_2
t_1 = max([t_A, t_B, t_C_prime])
t_out = t_1+t_2+t_3+2*N_ABC

N_ref = N_ABC

coal_ABC = N_ref/N_ABC
coal_AB = N_ref/N_AB
t_upper = t_3-cutpoints_ABC(n_int_ABC, 1/N_ABC)[-2]
t_AB = t_2/N_ref

cut_AB = t_1+cutpoints_AB(n_int_AB, t_AB, coal_AB)*N_ref
cut_ABC = t_1+t_2+cutpoints_ABC(n_int_ABC, coal_ABC)*N_ref

obs_dct = get_obs_state_dct()
dct = {obs_dct[i]: i for i in range(len(obs_dct))}

####################### Add demography #######################

print('Simulating demography')

demography = msprime.Demography()
demography.add_population(name="A", initial_size=N_AB, default_sampling_time=t_1-t_A)
demography.add_population(name="B", initial_size=N_AB, default_sampling_time=t_1-t_B)
demography.add_population(name="B_anc", initial_size=N_AB, initially_active=False)
demography.add_population(name="C", initial_size=N_BC, default_sampling_time=t_1+t_2-t_C)
demography.add_population(name="D", initial_size=N_AB, default_sampling_time=t_1-t_1)
demography.add_population(name="AB", initial_size=N_AB)
demography.add_population(name="ABC", initial_size=N_ABC)
demography.add_population(name="ABCD", initial_size=N_ABC)
demography.add_admixture(time = t_1-t_m, derived="B", ancestral=["B_anc", "C"], proportions=(1-m, m))
demography.add_population_split(time=t_1, derived=["A", "B_anc"], ancestral="AB")
demography.add_population_split(time=t_1+t_2, derived=["AB", "C"], ancestral="ABC")
demography.add_population_split(time=t_1+t_2+t_3, derived=["ABC", "D"], ancestral="ABCD")

ts = msprime.sim_ancestry(
    {"A": 1, "B": 1, "C": 1, 
     "D": 1
    }, 
    demography=demography, 
    recombination_rate=r,
    sequence_length=n_sites,
    ploidy=1, 
    random_seed=seed,
    num_replicates=10
)

####################### Add mutations #######################

print('Adding mutations')

E = []

for ts_individual in ts:

    mutated_ts = msprime.sim_mutations(ts_individual, rate=mu, random_seed=seed)

    nochange_lst = [dct['AAAA'], dct['CCCC'], dct['TTTT'], dct['GGGG']]
    sim_genome = np.random.choice(nochange_lst, n_sites)

    mut_lst = []
    mut_loc = []
    for variant in mutated_ts.variants():
        mut_loc.append(variant.site.position)
        mut_lst.append(''.join([variant.alleles[i] for i in variant.genotypes]))

    for i in range(len(mut_loc)):
        sim_genome[int(mut_loc[i])] = dct[mut_lst[i]]

    E.append(sim_genome)

####################### Optimization #######################

print('Optimizing')

t_A = t_A*mu
t_B = t_B*mu
t_C = t_C*mu
t_2 = t_2*mu
t_upper = t_upper*mu
t_out = t_out*mu
t_m = t_m*mu
N_AB = N_AB*mu
N_BC = N_BC*mu
N_ABC = N_ABC*mu
r = r/mu

transitions, emissions, starting, hidden_states, observed_states = trans_emiss_calc_introgression(
    t_A, t_B, t_C, t_2, t_upper, t_out, t_m,
    N_AB, N_BC, N_ABC,
    r, m, n_int_AB, n_int_ABC)

loglik = loglik_wrapper_par(transitions, emissions, starting, E)

write_list([-1, t_A, t_B, t_C, t_2, t_upper, t_m, N_AB, N_BC, N_ABC, r, m, loglik, 0], '../results/sim_{}_{}_{}_{}.csv'.format(n_int_AB, n_int_ABC, seed, model))

fct_big = 25
fct_sma = 2
t_sd_A = t_A/10000
t_sd_B = t_B/fct_sma
t_sd_C = t_C/fct_big
t_sd_2 = t_2/fct_big
t_sd_upper = t_2/fct_sma
t_sd_m = t_m/fct_sma
N_sd_AB = N_AB/fct_sma
N_sd_BC = N_BC/fct_sma
N_sd_ABC = N_ABC/fct_sma
r_sd = r/fct_sma

np.random.seed(seed)
t_init_A = np.random.uniform(t_A-t_sd_A, t_A+t_sd_A)
t_init_B = np.random.uniform(t_B-t_sd_B, t_B+t_sd_B)
t_init_C = np.random.uniform(t_C-t_sd_C, t_C+t_sd_C)
t_init_2 = np.random.uniform(t_2-t_sd_2, t_2+t_sd_2)
t_init_upper = np.random.uniform(t_upper-t_sd_upper, t_upper+t_sd_upper)
t_init_m = np.random.uniform(0, min([t_B-t_sd_B, (t_C-t_sd_C)-(t_2+t_sd_2)]))
N_init_AB = np.random.uniform(N_AB-N_sd_AB, N_AB+N_sd_AB)
N_init_BC = np.random.uniform(N_BC-N_sd_BC, N_BC+N_sd_BC)
N_init_ABC = np.random.uniform(N_ABC-N_sd_ABC, N_ABC+N_sd_ABC)
r_init = np.random.uniform(r-r_sd, r+r_sd)
m_init = np.random.uniform(0.0001, 0.5)

dct = {
    't_A':     [t_init_A,     t_A-t_sd_A, t_A+t_sd_A], 
    't_B':     [t_init_B,     t_B-t_sd_B, t_B+t_sd_B], 
    't_C':     [t_init_C,     t_C-t_sd_C, t_C+t_sd_C], 
    't_2':     [t_init_2,     t_2-t_sd_2, t_2+t_sd_2], 
    't_upper': [t_init_upper, t_upper-t_sd_upper, t_upper+t_sd_upper], 
    't_m':     [t_init_m,     0, min([t_B-t_sd_B, (t_C-t_sd_C)-(t_2+t_sd_2)])], 
    'N_AB':    [N_init_AB,    N_AB-N_sd_AB,  N_AB+N_sd_AB], 
    'N_BC':    [N_init_BC,    N_BC-N_sd_BC,  N_BC+N_sd_BC], 
    'N_ABC':   [N_init_ABC,   N_ABC-N_sd_ABC,  N_ABC+N_sd_ABC], 
    'r':       [r_init,       r-r_sd,  r+r_sd],
    'm':       [m_init,       0.0001,  0.5]
    }

x = ""
for i in dct:
    exec("x = %s" % (i))
    print(i, x > dct[i][1], x < dct[i][2])
    if (not (x > dct[i][1])) or (not (x < dct[i][2])):
        sys.exit()

print(dct)

dct2 = {'n_int_AB':n_int_AB, 'n_int_ABC':n_int_ABC}
res = optimizer_introgression(
    optim_params = dct, 
    fixed_params = dct2, 
    V_lst = E, 
    res_name = f'../results/sim_{n_int_AB}_{n_int_ABC}_{seed}_{model}.csv', 
    header = False,
    method = algorithm
    )


