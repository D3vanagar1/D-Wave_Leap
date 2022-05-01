# Followed Dwave Ocean Boolean NOT Gate example
# https://docs.ocean.dwavesys.com/en/latest/examples/not.html#not

# install necessary libraries
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
sampler = EmbeddingComposite(DWaveSampler())

'''
For a boolean not gate we must formluate it as a BQM
    * given M variables each var can have a value of 0 or 1, systems attempts to find values of var to minimize sum(q_i * x_i) + sum(q_ij * x_i * x_j), i<j
    * q_i and q_ij are linear and quadratic coefficents
Then represent it with a penalty function:
    2xz - x - z + 1 (where z = not(x))

Formulate Problem as QUBO (https://docs.ocean.dwavesys.com/en/latest/concepts/index.html#term-QUBO)
    * Note: QUBO matrix Q must be upper-diagonal (this happens anyways as we can always x2x1 = x1x2)
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

## Minor-Embedding a NOT Gate
'''
Minor Embedding maps x and z to the indexed qubits of the D-Wave QPU
Will need to use the FixedEmbeddingComposite library instead of EmbeddingComposite
in order to select a specific node of the sampler and look at its coupled qubits

* Mapping the NOT problem (2 linear coefficents and single quadratic coefficent)
  to biases on the D-Wave 2000Q's qubits 0 and 4 and coupling
  * NOT gate minor embedded into topmost left unit cell of D-Wave 2000Q QPU
  * x1 and x2 are minor embedded as qubits 0 and 4
  * q1,q2 = -1,-1 are the biases of each qubit and q12 = 2 is the coupling strength
  * images/NOT_gate/Embedding_Chimera_NOT.png
'''

sampler_embedded = FixedEmbeddingComposite(sampler, {'x': [0], 'z': [4]})
print(sampler_embedded.adjacency["x"])
sample_set_embedding = sampler_embedded.sample_qubo(Q, num_reads=7000, label='Sample - NOT Gate')
print(sample_set_embedding)

