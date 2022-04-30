
# install necessary libraries
from dwave.system import DWaveSampler, EmbeddingComposite
sampler = EmbeddingComposite(DWaveSampler())

'''
For a boolean not gate we must formluate it as a BQM
    * given M variables each var can have a value of 0 or 1, systems attempts to find values of var to minimize sum(q_i * x_i) + sum(q_ij * x_i * x_j), i<j
    * q_i and q_ij are linear and quadratic coefficents
Then represent it with a penalty function:
    2xz - x - z + 1 (where z = not(x))

Formulate Problem as QUBO (https://docs.ocean.dwavesys.com/en/latest/concepts/index.html#term-QUBO)
    * to map penalty function to BQN must reformulate our function
        - can drop -1 (const) as the the function's values are just shifted by -1
        - reorder in standard QUBO format:
            -x1 - x2 + 2x1x2 (z=x2, x=x1)
            q1 = q2 = -1, q12 = 2
            ==> Q = [[-1, 2],
                    [0, -1]]

Solve problem with sampling
'''

# QUBO coefficent matrix
Q = {('x', 'x'): -1, ('x', 'z'): 2, ('z', 'x'): 0, ('z', 'z'): -1}
sample_set = sampler.sample_qubo(Q, num_reads=5000, label='Samples - NOT Gate')
print(sample_set)

# check if correct
print(sample_set.first)
