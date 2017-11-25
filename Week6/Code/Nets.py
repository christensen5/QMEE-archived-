"""
Plot a network (using networkx) of the provided QMEE collaboration data.
"""

# IMPORTS
import networkx as nx
import scipy as sc
import matplotlib.pyplot as plt
import csv

# IMPORT PROVIDED NETWORK DATA FROM FILES
nodes_file = open("../Data/QMEE_Net_Mat_nodes.csv", 'rb')
edges_file = open("../Data/QMEE_Net_Mat_edges.csv", 'rb')

nodes = []
csvread = csv.reader(nodes_file)
csvread.next()  # skip header row
for row in csvread:
    nodes.append(tuple(row))
nodes_file.close()

edges = sc.array
csvread = csv.reader(edges_file)
AdjNames = list(csvread.next())
tmp = csvread.next()
tmp = map(int, tmp)
Adj = sc.asmatrix(tmp)
for row in csvread:
    tmp = map(int, row)
    Adj = sc.append(Adj, [tmp], axis = 0)
edges_file.close()


###### PLOTTING #######
plt.close('all')

G = nx.Graph(Adj)
nx.draw_circular(G)
plt.show()

