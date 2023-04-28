
from scipy import stats

import argparse

from write_tree import *
from TREE import *
from get_tree import *
from mcmc_movements import *
import time
from multiprocessing import Queue
from multiprocessing import Process
from multiprocessing import SimpleQueue
from multiprocessing import Pool
def task(tree, searches, m_h, pi, p_lamda, ifFixed, burnin):
  trees = MCMC(searches, tree, m_h, pi, p_lamda, ifFixed, burnin)
  return trees

def simulate_main(treeFile, snvTable, overlapped, reveal, result, searches, initAlpha = 0.01, initBeta = 0.2, initSigma = 0.05, restart = 10, ifFixed = 0, burnin = 2000, SNVlist = None, SNVcells = None):

  imputedFP =initAlpha
  imputedFN = initBeta
  imputedMiss = -1  

  treeText = result + ".newick"
  f_text = open(treeText, "w")
  f_log = open(result + ".log", "w")
  eval = result + ".Prob"
  f_eval = open(eval, "w")

  g_file = result + ".G"
  f_g = open(g_file, "w")
  snv_file = result + ".snv"
  f_snv = open(snv_file, "w")
  f_cell = open(result + ".cell", "w")
  topTree = None
  initTree, imputedMiss = readNewSimulatedTree(treeFile, snvTable, overlapped, reveal, initAlpha, initBeta, initSigma, SNVlist, SNVcells)

  queue = Queue()
  processes = []
  topTrees = []
  pi = 0.1#error rate move
  p_lamda = 0.1  # sigma move
  m_h = 0.5 #HM moves 

  MCMCargs = [(deepcopy(initTree), searches, m_h, pi, p_lamda, ifFixed, burnin) for _ in range(restart)]
  with Pool(processes=restart) as pool:
      topTrees = pool.starmap(task, MCMCargs)
  topTree = topTrees[0]
  for i in range(restart):
    
    f_log.write("restart: " + str(i) + "\n")
    tree = topTrees[i]
    if(tree['p'][0] > topTree['p'][0]):
      topTree = tree
    f_log.write("Iteration " + str(tree['id']) + "\n")
    f_log.write("Current probability is %s\n" % tree['p'][0])
    f_log.write("Restart\tTree_ID\tinitialFP\tinitalFN\tinitGamma\tinitialSigma\t" +
    "imputedFP\timputedFN\timputedMiss\testimatedFP\t" +"estimatedFN\testimatedMiss\testimatedSigma\t" +
    "Prob\t" + 
    "P_DG\tP_MG\tP_MC\tP_alpha\tP_beta\tP_gamma\tP_sigma\tP_G\n")
    f_log.write("{}\t{}\t{}\t{}\t{}\t{}".format(str(i), str(tree['id']), str(initAlpha), str(initBeta), str(imputedMiss), str(initSigma)))
    f_log.write("\t{}\t{}\t{}\t".format(str(imputedFP), str(imputedFN),str(imputedMiss) ))
    f_log.write("{}\t{}\t{}\t".format(str(tree['theta'][0]), str(tree['theta'][1]), str(tree['theta'][2])))
    f_log.write(str(tree['sigma']) + "\t")
    f_log.write("\t".join([str(p_) for p_ in tree['p']]) + "\n")
    writeTree(f_log, tree)
    f_log.write(tree['treeText'])

  G = topTree['G']
  for row in G:
      f_g.write(' '.join([str(a) for a in row]) + "\n")  
  f_text.write(topTree['treeText'])
  for s in topTree['snv']:
    f_snv.write(str(s) + "\t" + str(topTree['snv'][s]) + "\t" +str(topTree['snv_loss'][s]) +"\n")
  f_eval.write(str(topTree['p'][0]))
  for c in range(len(topTree['cells'])):
    f_cell.write(str(c) + "\t" + str(topTree['cells'][c]) + "\n")
             
  return 0


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-tree', help = "Required: tree structure", type = str)
  parser.add_argument('-D', help = "Required:  SNV table", type = str)
  parser.add_argument('-overlap', help = "Required:  file contain overlapped SNV and CNA", type = str, default = "None")
  parser.add_argument('-out', help = "Required: Output file", type = str)
  parser.add_argument('-alpha', help = "False positive value", type = float, default= 0.01)
  parser.add_argument('-beta', help = "False negative value", type = float, default = 0.2)
  parser.add_argument('-sigma', help = "Standar deviation for snv proportion and cn proportion", type = float, default = 0.05)
  parser.add_argument("-itr", help="Number of iteraton", type=int, default=100000)
  parser.add_argument("-reveal", help = "Reveal file: where SNV should be placed", type  = str, default="None")
  parser.add_argument("-restart", type = int, default=5)
  parser.add_argument("-fixCell", help="Whether the cell placement is fixed", type = int, default = 0)
  parser.add_argument("-burnin", help = "Burn it", default=2000, type=int)
  parser.add_argument("-SNVcell", help="SNV cells name and id file", type=str, default=None)
  parser.add_argument("-SNVid", help="SNV id and name", type=str, default=None)
  args = parser.parse_args()

  start_time = time.time()
  simulate_main(args.tree, args.D, args.overlap, args.reveal, args.out, args.itr, args.alpha, args.beta, args.sigma,  args.restart, args.fixCell, args.burnin, args.SNVid, args.SNVcell)
  print("--- %s seconds ---" % (time.time() - start_time))
