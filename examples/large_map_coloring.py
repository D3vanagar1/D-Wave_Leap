# Following steps from https://docs.ocean.dwavesys.com/en/latest/examples/map_kerberos.html#map-kerberos
# solve CSP (constraint satisfaction problem) to demo Ocean's hybrid sampler (KerberosSampler)


'''
* Formulate Problem
    * create text file (./text/usa.adj) that contains all US states and their adjacencies
    * convert file into a graph with networkx library
* Solve with Sampling
    *
'''
import networkx as nx
import dwave_networkx as dnx
from hybrid.reference.kerberos import KerberosSampler
import matplotlib.pyplot as plt
# represents states as vertices, neighbors share edges
G = nx.read_adjlist('./text/usa.adj', delimiter=',')

colouring = dnx.min_vertex_color(G, sampler=KerberosSampler(), chromatic_ub=4, max_iter=10, convergence=3)
set(colouring.values())

node_colours = [colouring.get(node) for node in G.nodes()]
if dnx.is_vertex_coloring(G, colouring):
    #nx.draw(G,pos=nx.shell_layout(G, nlist=[list(G.nodes)[x:x+10] for x in range(0,50,10)] + [[list(G.nodes)[50]]]), with_labels=True, node_color=node_colours, node_size=400, cmap=plt.cm.rainbow)
    # doesn't work "list index out of range"
plt.show()

