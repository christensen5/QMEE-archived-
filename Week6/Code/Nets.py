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
Adj = sc.asmatrix(csvread.next())
for row in csvread:
    Adj = sc.append(Adj, [row] , axis = 0)
edges_file.close()


###### PLOTTING #######
plt.close('all')

pos = nx.circular_layout(AdjNames)

G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(Adj)

