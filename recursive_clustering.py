#!/usr/bin/env python
# corrected index 0
# recursive script to cluster networks until quality measure is lower than cutoff
# written by Jose Flores-Canales, 08/2016
# used for DREAM challenge 11

import os
import sys
import atexit 
import subprocess
import networkx as nx
from collections import Counter
import numpy as np
import pickle
import re
import time
import glob
#import multiprocessing
try:
  from pylab import *
  import matplotlib as mpl
  #mpl.use('Agg')
  import matplotlib.pyplot as plt
except:
  print "Problems importing Matplotlib and Pylab modules"

USAGE="""
recursive_clustering.py edgefile gmin
"""

proc_list = []

###################################################################################
# call_loop inquires jobs status in the cluster
###################################################################################
def call_loop(jobids, pwd_dir):
  flags = {}
  n = 1
  while True:
    #print "Loop: ", n
    for idx, pathid in jobids.iteritems():
      cmd = ['qstat', '-a'] 
      proc = subprocess.Popen(cmd, shell = False, stdout=subprocess.PIPE)
      proc_list.append(proc)
      (out, err) = proc.communicate()
      if not re.search(idx, out) and not (idx in flags):
        flags[idx] = True
        print "Parent path_dir :", pwd_dir, "Child pathid: ", pathid
        os.chdir('%s' % (pathid))
        edgelist = find_edgelist(pathid)
        if edgelist:
          cmd = ['python', 'recursive_clustering.py', edgelist, 'gmin.txt', pathid]
          proc1 = subprocess.Popen(cmd, shell = False, stdout=subprocess.PIPE)
          proc_list.append(proc1)
          (out1, err1) = proc1.communicate()
    if len(jobids) == len(flags):
      print " Cleaning... "
      break
    time.sleep(30)
    n += 1

def write_log(path, buffer):
  filename = path + 'clustering.log'
  if os.path.exists(filename):
    with open(filename, 'a') as f:
      for bitem in buffer:
        for i in bitem: 
          print >> f, i 
        print >> f, '\n'
  else:
    with open(filename, 'w') as f:
      for bitem in buffer:
        for i in bitem: 
          print >> f, i 
        print >> f, '\n'

def write_dir_info(path, n_nodes):
  filename = path + '/info.folder.txt'
  f = open(filename, 'w')
  f.write('%d' % n_nodes)
  f.close() 

def read_dir_info(path):
  filename = path + '/info.folder.txt'
  with open(filename, 'r') as f:
    lines = f.readlines()
  return int(lines[0].strip().split()[0])

def find_edgelist(path):
  files_e = glob.glob(path + '/' + 'edgelist*.txt')
  if len(files_e) != 0:
    return files_e[0]
  else:
    return '' 

def Sum_intercomm_weights(G, community, i, j):
  S = 0
  if i == j:
    for ii in xrange(0, len(community[i])):
      for jj in xrange(ii, len(community[i])):
        n, m = community[i][ii], community[i][jj]
        if m in G[n]:
          S += G[n][m]['weight'] 
  else:
    for n in community[i]:
      for m in community[j]:
        if m in G[n]:
          S += G[n][m]['weight'] 
  return S

###################################################################################
# prints PBS scripts, modify accordingly to your system, queue is specified at
# sumbission 
###################################################################################
def print_runfile(k, dirname):
  run_lines =['#!/bin/sh\n',
  '### Job name\n',
  '#PBS -N ',  # add name
  '### Declare job non-rerunable\n',
  '#PBS -r n\n',
  '### Output files\n',
  '#PBS -j oe\n',
  '### Mail to user\n',
  '#PBS -m ae \n',
  '### Queue name (n1, n2, n4, n16, n60, g2, g4)\n',
  '#PBS -q n24\n',
  '### Walltime limit (hh:mm:ss)\n',
  '#PBS -l walltime=360:00:00\n',
  'echo Working directory is $PBS_O_WORKDIR\n',
  'cd $PBS_O_WORKDIR\n',
  'echo Running on host `hostname`\n',
  'echo Time is `date`\n',
  'echo Directory is `pwd`\n',
  'echo This job runs on the following processors:\n',
  'echo `cat $PBS_NODEFILE`\n',
  'NPROCS=`wc -l < $PBS_NODEFILE`\n',
  'mpirun -np $NPROCS  ./comm_csa.mpi.weighted.x\n']
  f = open(dirname+'/'+'run.sh','w')
  i = 1
  comm = 'comm'+str(k)
  for line in run_lines:
    if i == 4:
      line_n = '%s\n' % comm 
      f.write(line_n)
    f.write(line)
    i += 1
  f.close()

def print_csafile(k, atom, dirname, edgelistfile):
  csa_lines = ['50 100           nconf,nbankm\n',
  '      jstart,jend\n', #add i j 
  '999999 999999    timem,iterm\n',
  '5 10 5 5 2 ', 
  '      n1,n2,n3,n4,is1,is2 : n1 -> br, n2-> bb, n3->newconf_residue, n4->mutation\n', # add number of 0.4*nodes
  '5 3 5            nran0,nran1,irr\n',
  '50 10            nseed,nadd\n',
  '70000 2.0 5.0    ntotal,cut1,cut2\n',
  '-999999.99       estop for N=55\n',
  '1 2 0 0          icmax,iucut,irestart,ifbank\n',
  '0 1              min_type dist_type(0->variation_of_info, 1->mutual_info, 2-> common_pair)\n']
  if atom*0.4 < 44: 
    g_atom = 44
  else:
    g_atom = int(np.floor(atom*0.4))
  
  f = open(dirname+'/'+'csa.in','w')
  i = 1
  k = 1 # this can be modified to run multiple modcsa runs 
  for line in csa_lines:
    if i == 2:
      line_n = '%d %d' % (k,k)
      f.write(line_n)
    elif i == 5:
      line_n = '%d' % g_atom
      f.write(line_n)
    f.write(line)
    i += 1
  line = '%s\n' % edgelistfile
  f.write(line)
  f.close()

###################################################################################
# main function of this script, requires edgelist, gmin, and current working
# directory path 
###################################################################################
def main_run(edgelist, gmin):
  buffer = []

  with open(edgelist,'r') as networkFile:
    lines_data = networkFile.readlines() 

  Gnew = nx.Graph(name = 'Network')
 
  for line in lines_data :
    if not line.startswith(('\n','\t','#')):
      try:
        g1,g2, w = line.split('\t')[:3]
        Gnew.add_edge(int(g1), int(g2), weight = float(w))
      except:
        continue
 
  ########### PRINT INFO ##############################################
 
  # load membership output from modcsa 
  with open(gmin,'r') as f:
    lines = f.readlines() 
  values = []
  node = 1 
  for line in lines:
    if not line.startswith(('#')):
      nums = line.strip().split()
      values.append(int(nums[0]))
      Gnew.node[node]['membership'] = int(nums[0])
      node += 1
    else:
      nums = line.strip().split()
      Q = float(nums[3]) 

  # count the number of nodes per community
  # generate dictionary [community_i] of sizes of communities
  counted_membership =  dict(Counter(values))
  buffer.append(["Counted Modules: d[i] = size_i: ", counted_membership])

  # generate dictionary [community_i] of lists of node indexes
  communities = {}
  node = 1 
  for value in values:
    if value in communities.keys():
      communities[value].append(node)
    else:
      communities[value] = [node]
    node += 1

  # get the size of communities in the network
  N_comm = len(list(communities))
  buffer.append(["Number of Modules: ", N_comm])

  # generate dictionary of labels
  labels = {}
  for v in list(communities):
    labels[v-1] = str(v)
  #  print v, communities[v]

  # generate matrix of interactions between communities
  matrix = np.zeros((N_comm,N_comm),dtype = np.float)
  for i in xrange(1,N_comm+1):
    for j in xrange(i, N_comm+1):
      S = 0
      matrix[i-1,j-1] = Sum_intercomm_weights(Gnew, communities, i, j)
      matrix[j-1,i-1] = matrix[i-1,j-1] 

  ########################################################################################################
  # print general information about each node (community)
  nodes_weight, nodes_len = [], []
  for i in xrange(N_comm):
    nodes_weight.append(matrix[i,i])
    nodes_len.append(len(communities[i+1]))
  buffer.append(["Nodes Weight:", nodes_weight, "Number of nodes", len(nodes_weight)])
  buffer.append(["Nodes Size:", nodes_len, "number of nodes", len(nodes_len)])

  # reorder the matrix of interactions according to weight and size of nodes:
  axes_weight = np.argsort(np.array(nodes_weight),axis=-1)[::-1]
  axes_len = np.argsort(np.array(nodes_len),axis=-1)[::-1]
  buffer.append(["Axes weight: ", axes_weight])
  buffer.append(["Axes len: ", axes_len])

  community_max_weight = int(labels[axes_weight[0]])
  community_max_len = int(labels[axes_len[0]])

  #print "old labels", labels.values()
  labels_sorted_by_len = np.array(labels.values())[axes_len]
  labels_community = []
  for i in xrange(N_comm):
  #  labels[i] = labels_sorted_by_len[i]
    labels_community.append(int(labels_sorted_by_len[i]))
  #print "new labels", labels.values()

  buffer.append(["Community with the largest weight", community_max_weight])
  buffer.append(["Community with the largest number of nodes", community_max_len])

  nodes_weight_sorted = np.array(nodes_weight)[axes_weight]
  nodes_len_sorted = np.array(nodes_len)[axes_len]
  buffer.append(["Nodes weight sorted", nodes_weight_sorted])
  buffer.append(["Nodes len sorted", nodes_len_sorted])

  matrix_sorted_by_weight = matrix[axes_weight,:][:,axes_weight]
  matrix_sorted_by_len = matrix[axes_len,:][:,axes_len]

  #### PRINT matrices of intra- and inter-sum of weights of communities

  np.savetxt('matrix_sorted_by_weight.txt',matrix_sorted_by_weight,fmt='%10.1f', delimiter='',newline='\n')
  np.savetxt('matrix_sorted_by_len.txt',matrix_sorted_by_len,fmt='%10.1f', delimiter='',newline='\n')

  ###################################################################################
  # Q > -0.35, not modular enough  
  ###################################################################################
  if Q > -0.001:
    return 
  buffer.append(["Q:", Q])
  ###################################################################################
  # find number of cores per node 
  # queues are n8, n16, or n24
  ###################################################################################
  #n_cores = multiprocessing.cpu_count()
  n_cores = 24 

  ###################################################################################
  # create a directory for each community
  # define the constant to remove from the edge weights
  ###################################################################################

  K = 0.0
  pwd_dir = os.getcwd()
  pwd_dir += '/'
  write_log(pwd_dir, buffer) 

  path_exe = '/home/jfloresc/dream11/subchallenge1/map_module.py' 

  write_dir_info(pwd_dir, len(Gnew))
  list_jobids = {} 
  for i in xrange(0,N_comm):
    print "printing community:", labels_community[i]
    nodes_max = [gene1 for (gene1, w) in Gnew.nodes(data=True) if w['membership'] == labels_community[i]]
    Gtemp = Gnew.subgraph(nodes_max)
    #nameout = 'community'+str(i+1)+'.gn' 
    #nx.write_gpickle(Gtemp,nameout,pickle.HIGHEST_PROTOCOL)
    dirname = 'comm' + str(i+1)
    if not os.path.isdir(dirname):
      os.makedirs(dirname)
    
    nameout = 'community'+str(i+1)+'.txt'
    namemapped = 'edgelist'+str(i+1)+'.txt'
    path_nameout = dirname + '/' + nameout
    path_dirname = pwd_dir + dirname 

    nx.write_weighted_edgelist(Gtemp, path_nameout, delimiter='\t')
    write_dir_info(path_dirname, len(Gtemp))
    
    if len(Gtemp) > 100:
      os.system('cp /home/jfloresc/modcsa/src/comm_csa.mpi.weighted.x recursive_clustering.py %s' % (dirname))
      os.system('python %s %s %s %s %10.6f'% (path_exe,dirname, nameout, namemapped, K))
      print_csafile(i+1, len(Gtemp), dirname, namemapped) 
      print_runfile(i+1, dirname) 

      os.chdir('%s' % (path_dirname))
      ###################################################################################
      # queues are n8, n16, or n24
      ###################################################################################
      #cmd_t1 = '-q n%d' % (n_cores)
      cmd_t2 = '%s' % ('run.sh')
      cmd = ['qsub', cmd_t2]
      proc = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE)
      proc_list.append(proc)
      (out, err) = proc.communicate()
      jobid = out.split('.')[0]
      list_jobids[jobid] = path_dirname 
      os.chdir('%s' % (pwd_dir))

  call_loop(list_jobids,pwd_dir)

###################################################################################
# attempts to kill all processes  
###################################################################################
def cleanup():
  for p in proc_list:
    try:
      p.kill()
    except:
      pass

def main():
  if len(sys.argv)<2:
    print USAGE
    sys.exit()
  edgelist = sys.argv[1]
  gmin = sys.argv[2]
  main_run(edgelist, gmin)
  atexit.register(cleanup)

if __name__ == '__main__':
    main()

