# Following Dwave Ocean Boolean AND Gate example
# https://docs.ocean.dwavesys.com/en/latest/examples/and.html#and
# For embedding need to know about Chimera: https://docs.ocean.dwavesys.com/en/latest/concepts/topology.html

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
import dwave.inspector
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
sampler = DWaveSampler()
sampler_embedded = EmbeddingComposite(sampler)

# QUBO matrix
Q = {('x1', 'x2'): 1, ('x1', 'z'): -2, ('x2', 'z'): -2, ('z', 'z'): 3}
sample_set = sampler_embedded.sample_qubo(Q, num_reads = 10000, label = 'Samples - AND Gate')
print(sample_set)
# Note somtimes some columns have the same values everywhere except for rightmost column (chain_breaks)


## Non-automated Minor-Embedding
'''
For the NOT gate we had to minor-embed a K2 graph (straight line)
But for an AND gate we have to embed a K3 graph (triangle) on D-Wave 2000Q system. Will require chaining qubits
(Note: wouldn't need chaining on Pegasus system as it can minor-embed K3)
* Minor-Embedding of K3 graph (AND gate)
    * looking at Chimera topology (images/AND_gate/Chimera_Unit_Cell.png) cannot connect 3 qubits in closed loop
    * can make closed loop of 4 qubits (ex qubits 0,3,4,7)
    * create a chain of 2 qubits to represent a single variable => convert square to triangle
        * (process found in images/AND_gate/Embedding_Chimera_AND.png)
        * strength of 0-4 coupler must be set to strongly correlate the 2 qubits so that most solutions have single z value
            * Chain strength (https://docs.ocean.dwavesys.com/en/latest/concepts/embedding.html#concepts-chain-strength)
            * likely explains chain breaks seen earlier. Qubits in a chain with different values
        * use FixedEmbeddingComposite library for manual minor-embedding
'''

embedding = {'x1': {3}, 'x2': {7}, 'z': {0,4}}
sampler_minor_embedded = FixedEmbeddingComposite(sampler, embedding)
print(sampler_minor_embedded.adjacency)
sample_set_embedding = sampler_minor_embedded.sample_qubo(Q, num_reads = 10000, label = 'Samples - AND Gate')
print(sample_set_embedding)

# view solution on QPU
# graphical view of og BQM representation and its embedded representation (images/AND_gate/Inspector_AND_gate.png)
dwave.inspector.show(sample_set_embedding)
