# Following Dwave Ocean Boolean AND Gate example
# https://docs.ocean.dwavesys.com/en/latest/examples/and.html#and

'''
* Formulate AND Gate as BQM
    * Form objective function
        E = sum(h_i * s_i) + sum(J_ij * s_i * s_j) , i<j
        h_i are biases , J_ij couplings between spins
    * CS equivalent is the QUBO
* Create Penalty Func for AND
    * x1 ^ x2 <=> z (x1 , x2 are inputs)
    FUNC: f == x1x2 - 2(x1 + x2)z + 3z
    CHECK:
       z | x2 | x1 || f | T/F
       0 | 0  | 0  || 0 | T
       0 | 0  | 1  || 0 | T
       0 | 1  | 0  || 0 | T
       0 | 1  | 1  || 1 | F
       1 | 0  | 0  || 3 | F
       1 | 0  | 1  || 1 | F
       1 | 1  | 0  || 1 | F
       1 | 1  | 1  || 0 | T
    all cases when f = 0 its correct and non zero values of f are false => function works
* Formulate as QUBO problem
    * E = 3 * x3 + x1x2 - 2 * x1x3 -2 * x2x3
    * z = x3, q1 = 3, q12 = 1, q13 = -2, q23 = -2
    ==> Q = [[0,1,-2],
             [0,0,-2],
             [0,0,3]]
'''

#import necessary libraries
from dwave.system import DWaveSampler, EmbeddingComposite
sampler = DWaveSampler()
sampler_embedded = EmbeddingComposite(sampler)

# QUBO matrix
Q = {('x1', 'x2'): 1, ('x1', 'z'): -2, ('x2', 'z'): -2, ('z', 'z'): 3}
sample_set = sampler_embedded.sample_qubo(Q, num_reads = 10000, label = 'Samples - AND Gate')
print(sample_set)
# Note somtimes some columns have the same values everywhere except for rightmost column (chain_breaks)
