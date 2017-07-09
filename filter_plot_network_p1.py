#!/usr/bin/env python
# written by Jose Flores-Canales, 08/2016
# used for DREAM challenge 11
 
import networkx as nx
import pickle
from pylab import *
 
with open('1_ppi_anonym_v2.txt','r') as networkFile:
    lines_data = networkFile.readlines() 

NetworkP3 = nx.Graph(name = 'Network P1')
 
for line in lines_data :
  if not line.startswith(('\n','\t','#')):
    g1,g2, w = line.split('\t')[:3]
    NetworkP3.add_edge(g1,g2,weight=float(w))
 
########### SAVE AND LOAD ##############################################
 
with open('1.gn', 'w') as output:
    pickle.dump(NetworkP3, output, pickle.HIGHEST_PROTOCOL)
 
NetworkP3 = pickle.load(open('1.gn', 'r'))
 
#### PRINT INFOS ON THE NETWORK ########################################
 
undirected_net = NetworkP3.to_undirected()
 
print 'General informations : \n', \
        nx.info(NetworkP3)
 
connected_components = nx.connected_components(undirected_net)
 
print 'Number of connected components : ', \
        len(list(connected_components) )

print 'Length of the connected components : ', \
       [len(c) for c in sorted(nx.connected_components(undirected_net), key=len, reverse=True)] 
 
#main_component = list(nx.connected_component_subgraphs(undirected_net))[0]
main_component = max(nx.connected_component_subgraphs(undirected_net), key=len)
 
print "Average shortest path in P3's main connected component : ", \
       nx.average_shortest_path_length(main_component)
 
print "Diameter of P3's main connected component : ", \
       nx.diameter(main_component)
 
NetworkMax = nx.Graph(data = main_component, name = 'Network Max P1')

with open('1max.gn', 'w') as output:
    pickle.dump(NetworkMax, output, pickle.HIGHEST_PROTOCOL)

nx.write_weighted_edgelist(NetworkMax, '1.maxcomponent.edgelist') 
 
NetworkMax = pickle.load(open('1max.gn', 'r'))
 
print 'General informations : \n', \
        nx.info(NetworkMax)
 
###########  D R A W  ##################################################
figure()
nx.draw(NetworkMax, node_size = 20)
show()
