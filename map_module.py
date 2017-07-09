#!/usr/bin/env python
# written by Jose Flores-Canales, 08/2016
# used for DREAM challenge 11

# modified mapping to begin the index from  index 1

import sys
import json


USAGE="""
map_module.py working_directory infile outfile K
working directory were input file is located
infile name of input file
outfile name of mapped edgelist file (written in working directory)
K parameter default 0
"""

def map_function(path, filename, nameout, K):
  with open(path+'/' + filename,'r') as f:
    lines = f.readlines()

  gene1_list, gene2_list, weights = [], [], []
  gene1_alt_list, gene2_alt_list = [], []

  for line in lines:
    nums = line.strip()
    nums = nums.split()
    gene1, gene2, w = int(nums[0]), int(nums[1]),float( nums[2])
    if w > K:  
      gene1_list.append(gene1)
      gene2_list.append(gene2)
      weights.append(w)
    else:
      gene1_alt_list.append(gene1)
      gene2_alt_list.append(gene2)

# number of edges in the edgelist file
  N = len(weights)
# get unique nodes in each list of genes  
  gene1_list_new, gene2_list_new = list(set(gene1_list)), list(set(gene2_list)) 
  n1, n2 = len(gene1_list_new), len(gene2_list_new) 

  gene1_t, gene2_t = sorted(gene1_list_new), sorted(gene2_list_new)
######
  gene_all = gene1_t + gene2_t
  gene_all_new = list(set(gene_all))
# get ordered list of unique nodes
  gene_all_new_s = sorted(gene_all_new)
###
  print "Number of genes1 :", n1
  print "Number of genes2 :", n2
  n_unique = len(gene_all_new_s)
  print "Number of unique genes :", n_unique 

  map_genes, inv_map_genes =  {}, {}
  c = 1 # changed to begin from 1 
  for i in gene_all_new_s:
    if i in map_genes.values():
      continue
    else:
      map_genes[c] = i

    if i in inv_map_genes.keys():
      continue
    else:
       inv_map_genes[i] = c

    c += 1

  with open(path+'/'+'map_genes.dat','w') as f:
    json.dump(map_genes,f)
    f.flush()
  with open(path+'/'+'inv_map_genes.dat','w') as f:
    json.dump(inv_map_genes,f)
    f.flush()


  outfile = open(path+'/'+nameout,'w')
  for i in xrange(N):
    gene1, gene2, w = gene1_list[i], gene2_list[i], float(weights[i]) 
    # modified to substract a constant K to all weighted edges
    line = "%d\t%d\t%9.6f\n" % (inv_map_genes[gene1],inv_map_genes[gene2],  w)
    outfile.write(line)
      
    #line = "%d\t%d\t%9.6f\t%d\t%d\n" % (inv_map_gene1[gene1], inv_map_gene2[gene2], w, gene1, gene2)
  line = "%d %d 0" % (n_unique,n_unique)
  outfile.write(line)
  outfile.close()


def main():
  if len(sys.argv)<4:
    print USAGE
    sys.exit()
  path = sys.argv[1]
  namein = sys.argv[2]
  nameout = sys.argv[3]
  K = float(sys.argv[4])
  map_function(path, namein, nameout, K)

if __name__ == '__main__':
    main()

