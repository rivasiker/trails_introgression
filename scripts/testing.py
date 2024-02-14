import numpy as np
from trails.cutpoints import cutpoints_ABC

mu = 2e-8
g = 25
N_AB = 10_000*2
N_ABC = 10_000*2
N_ref = N_ABC
t_1 = 60_000/g
t_2 = (600_000-t_1)/g
t_3 = (6_000_000-t_2)/g
t_A = t_1
t_B = t_1
t_C = t_1+t_2
r = 1e-8
seed = 1

n_int_ABC = 1

t_m = t_1-54_000/g
m = 0.05

t_upper = t_3-cutpoints_ABC(n_int_ABC, 1/N_ABC)[-2]
t_out = t_1+t_2+t_3+2*N_ABC

t_A = t_A*mu
t_B = t_B*mu
t_C = t_C*mu
t_2 = t_2*mu
t_upper = t_upper*mu
t_out = t_out*mu
t_m = t_m*mu
N_AB = N_AB*mu
N_ABC = N_ABC*mu
r = r/mu

fct_big = 25
fct_sma = 2
t_sd_A = t_A/fct_sma
t_sd_B = t_B/fct_sma
t_sd_C = t_C/fct_big
t_sd_2 = t_2/fct_big
t_sd_upper = t_2/fct_sma
t_sd_m = t_m/fct_sma
N_sd_AB = N_AB/fct_sma
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
    'N_ABC':   [N_init_ABC,   N_ABC-N_sd_ABC,  N_ABC+N_sd_ABC], 
    'r':       [r_init,       r-r_sd,  r+r_sd],
    'm':       [m_init,       0.0001,  0.5]
    }

for i in dct:
    exec("x = %s" % (i))
    print(i, x > dct[i][1], x < dct[i][2])

print(0, t_m, min([t_B-t_sd_B, (t_C-t_sd_C)-(t_2+t_sd_2)]))
print(t_m)