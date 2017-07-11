import pandas
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import pickle

ontology = pandas.read_table("070717.ontology.txt",sep="\t")

g = nx.DiGraph()
for i in ontology['Name'].tolist():
    g.add_node(i)
accession_name = dict(zip(ontology['Accession'].tolist(),ontology['Name'].tolist()))
nodes = []
for i in ontology.index:
    n1 = ontology.loc[i,'Accession']
    ns = ontology.loc[i,'Parents']
    if not pandas.isnull(ns):
        ns = ns.split(',')
        for j in ns:
            g.add_edge(accession_name[j],accession_name[n1])

pos=graphviz_layout(g,'dot')
nx.draw(g,pos=pos,node_size=7)
plt.rcParams["figure.figsize"] = (20,20)
plt.savefig('ontology.pdf')

g = g.edges()
data = {}
for i,j in g:
    if i not in data:
        data[i] = [j]
    else:
        data[i].append(j)

pickle.dump(data,open("070717.ontology.pickle","w"))
