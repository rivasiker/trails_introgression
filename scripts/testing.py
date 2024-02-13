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

t_m = t_1-40_000/g
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

np.random.seed(seed)
t_init_A = np.random.normal(t_A, t_A/5)
t_init_B = np.random.normal(t_B, t_B/5)
t_init_C = np.random.normal(t_C, t_C/5)
t_init_2 = np.random.normal(t_2, t_2/5)
t_init_upper = np.random.normal(t_upper, t_upper/5)
t_init_m = np.random.normal(t_m, t_m/5)
N_init_AB = np.random.normal(N_AB, N_AB/5)
N_init_ABC = np.random.normal(N_ABC, N_ABC/5)
r_init = np.random.normal(r, r/5)
m_init = m

# acc = 1
# while (t_init_B/2 <= t_init_m) or ((t_init_C/2-t_init_2*2) <= t_init_m):
#     np.random.seed(seed*10000+acc)
#     t_init_A = np.random.normal(t_A, t_A/5)
#     t_init_B = np.random.normal(t_B, t_B/5)
#     t_init_C = np.random.normal(t_C, t_C/5)
#     t_init_2 = np.random.normal(t_2, t_2/5)
#     t_init_upper = np.random.normal(t_upper, t_upper/5)
#     t_init_m = np.random.normal(t_m, t_m/5)
#     N_init_AB = np.random.normal(N_AB, N_AB/5)
#     N_init_ABC = np.random.normal(N_ABC, N_ABC/5)
#     r_init = np.random.normal(r, r/5)
#     m_init = m
#     acc += 1

# print(acc)

dct = {
    't_A':     [t_init_A,     t_init_A/2, t_init_A*2], 
    't_B':     [t_init_B,     t_init_B/2, t_init_B*2], 
    't_C':     [t_init_C,     t_init_C/2, t_init_C*2], 
    't_2':     [t_init_2,     t_init_2/2, t_init_2*2], 
    't_upper': [t_init_upper, t_init_upper/2, t_init_upper*2], 
    't_m':     [t_init_m,     0, min([t_init_B/2, t_init_C/2-t_init_2*2])], 
    'N_AB':    [N_init_AB,    N_init_AB/2,  N_init_AB*2], 
    'N_ABC':   [N_init_ABC,   N_init_ABC/2,  N_init_ABC*2], 
    'r':       [r_init,       r_init/5,  r_init*5],
    'm':       [m_init,       0.0001,  0.5]
    }

print(dct)