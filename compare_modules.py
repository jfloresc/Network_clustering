#!/usr/bin/env python
# written by Jose Flores-Canales, 08/2016
# used for DREAM challenge 11

#import mutual_information
import calc_NMI as NMI
import numpy as np
import os
import sys
import glob
import subprocess
import json
import matplotlib.pyplot as plt

USAGE="""
compare_modules.py Q2 10 Q0_K0.1
a) Prefix of directories: Q2, Q1, ...
b) N number of directories (suffix): Q2_1, Q2_2,...Q2_N 
c) subdirectory name
"""

min_len = lambda x,y: len(x) if len(x) < len(y) else len(y)

def byteify(input):
  if isinstance(input, dict):
    return {byteify(key):byteify(value) for key,value in input.iteritems()}
  elif isinstance(input, list):
    return [byteify(element) for element in input]
  elif isinstance(input, unicode):
    return input.encode('utf-8')
  else:
    return input

def read_q(path):
  namefile = path + '/q.txt'
  if os.path.isfile(namefile):
    with open(namefile, 'r') as f:
      data = json.load(f)
    data = byteify(data)
    return data
  else:
    print "ERROR: No modularity measure in modules"
    sys.exit()

def read_partition(pwd_dir):
  nodes, edges = [], []
  os.chdir('%s' % (pwd_dir))
  files = glob.glob('module_new*.txt')
  n_modules = len(files)
  if n_modules == 0:
    print "Error: No filtered modules at: ", os.getcwd()
    return []

  for i in xrange(n_modules):
    filename = files[i]
    with open(filename, 'r') as f:
      lines = f.readlines()

    nodes_i = set() 
    edges_i = []
    for line in lines :
      if not line.startswith(('\n','\t','#')):
        try:
          nums = line.strip().split()
          g1, g2, w = nums[0], nums[1], nums[3] 
          w = w.split('}')[0]
          edges_i.append([int(g1), int(g2), float(w)])
          nodes_i.add(int(g1))
          nodes_i.add(int(g2))
        except:
          continue
    nodes.append(nodes_i)
    edges.append(edges_i)

  return n_modules, nodes, edges 

def write_modules_stat(pwd, data, filename):
  namepng = pwd + '/' + filename + '.png' 
  xmin = 0
  xmax = 103
  bins = xrange(xmin, xmax + 1, 1)
  print "Max. Number of nodes in modules: ", max(data), "and Min.: ", min(data)
 #np.histogram(data, bins)
  plt.figure(1)
  plt.hist(data, bins)
  plt.title('Histogram of number of nodes per module')
  axes = plt.gca()
  axes.set_xlim([xmin, xmax])
  plt.savefig(namepng)
  plt.close()

def find_similar_modules(modules_i, modules_j):
  n_i, n_j = modules_i[0], modules_j[0] 
  module_pairs = []
  nodes_i, nodes_j = modules_i[2], modules_j[2]
  edges_i, edges_j = modules_i[3], modules_j[3]
  
  for i in xrange(n_i):
    for j in xrange(n_j):
      join = list(nodes_i[i] & nodes_j[j])
      if join:
        continue
      else:
         ratio = float(len(join)) / float(min_len(nodes_i[i], nodes_j[j]))
         if ratio > 0.8:
           module_pairs.append((i, j))

  return module_pairs       

def main_run(dir_prefix, dir_suffix_n, subdir):

  pwd = os.getcwd()
  modules_info = {} 
  for i in xrange(1, dir_suffix_n+1):
    dirname = dir_prefix + '_' + str(i)
    path_dir = pwd + '/' + dirname 
    path_sub = path_dir + '/' + subdir
    path_modules = path_sub + '/modules'
    path_modules_filtered = path_sub + '/t2'
      
    #print path_modules_filtered
    if not os.path.isdir(path_modules_filtered):
      print "Error: No Filtered Modules directory: t2", os.getcwd()
      sys.exit()

    n_modules, nodes, edges = read_partition(path_modules_filtered)  
    modules_info[i] = [n_modules, dirname, nodes, edges]
  
  pair_modules = {}
 
  m_list = [[] for i in xrange(dir_suffix_n)] 
  for i in xrange(1, dir_suffix_n ):
    for j in xrange(i + 1, dir_suffix_n + 1):
      #pairs_modules[(i, j)] = find_similar_modules(modules_info[i], modules_info[j])
      temp = find_similar_modules(modules_info[i], modules_info[j])
      m_list[i-1].append([j, temp])
      m_list[j-1].append([i, temp])

  filename = dir_prefix + '_' + subdir
  
  for i in xrange(1, dir_suffix_n ):
    print m_list[i-1]
  #for key, value in pairs_modules.iteritems():
    

def main():
  if len(sys.argv)<3:
    print USAGE
    sys.exit()
  dir_prefix = sys.argv[1]
  dir_suffix_n = int(sys.argv[2])
  subdir = sys.argv[3]
  main_run(dir_prefix, dir_suffix_n, subdir)

if __name__ == '__main__':
  main()
