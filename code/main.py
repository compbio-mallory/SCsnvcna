
from scipy import stats

import argparse

from write_tree import *
from data_structure import *
from TREE import *
from get_tree import *
from mcmc_movements import *
import time



def simulate_main(treeFile, snvTable, overlapped, reveal, result, searches, initAlpha = 0.001, initBeta = 0.2, initSigma = 0.05, restart = 10, ifFixed = 1):

  imputedFP =initAlpha
  imputedFN = initBeta
  imputedMiss = -1  

  treeText = result + ".treeText"
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
  initTree, imputedMiss = readNewSimulatedTree(treeFile, snvTable, overlapped, reveal, initAlpha, initBeta, initSigma)

  for i in range(restart):
    tree = deepcopy(initTree)
    #imputedFP, imputedFN = getError(tree.D, G) 
    
    pi = 0.1#error rate move
    p_lamda = 0.1  # sigma move
    #p_g = 0 # ground true data move
    m_h = 0.5 #HM moves 

    print("restart:", i)
    start_time = time.time()
    trees = MCMC(searches, tree, m_h, pi, p_lamda, ifFixed)
    print("--- %s seconds ---" % (time.time() - start_time))
    if not topTree:
      topTree = deepcopy(trees[0])
    else:
      if trees[0].p[0] > topTree.p[0]:
        topTree = deepcopy(trees[0])
    

    f_log.write("restart: " + str(i) + "\n")
    for tree in trees:
      f_log.write("Iteration " + str(tree.id) + "\n")
      f_log.write("Current probability is %s\n" % tree.p[0])
      f_log.write("Restart\tTree_ID\tinitialFP\tinitalFN\tinitGamma\tinitialSigma\t" +
      "imputedFP\timputedFN\timputedMiss\testimatedFP\t" +"estimatedFN\testimatedMiss\testimatedSigma\t" +
      "Prob\t" + 
      "P_DG\tP_MG\tP_MC\tP_alpha\tP_beta\tP_gamma\tP_sigma\tP_G\n")
      f_log.write("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(str(i), str(tree.id), str((tree.id)), str(initAlpha), str(initBeta), str(imputedMiss), str(initSigma)))
      f_log.write("\t{}\t{}\t{}\t".format(str(imputedFP), str(imputedFN),str(imputedMiss) ))
      f_log.write("{}\t{}\t{}\t".format(str(tree.theta[0]), str(tree.theta[1]), str(tree.theta[2])))
      f_log.write(str(tree.sigma) + "\t")
      f_log.write("\t".join([str(p_) for p_ in tree.p]) + "\n")
      writeTree(f_log, tree)
      f_log.write(genTreeText(tree))

  G = topTree.G
  for row in G:
      f_g.write(' '.join([str(a) for a in row]) + "\n")  
  f_text.write(genTreeText(topTree))
  for s in topTree.snvs:
    f_snv.write(str(s.id) + "\t" + str(s.edge) + "\t" +str(s.loss) +"\n")
  f_eval.write(str(topTree.p[0]))
  for c in range(len(topTree.cells)):
    f_cell.write(str(c) + "\t" + str(topTree.cells[c]) + "\n")
             
  return 0


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-tree', help = "Required: tree structure", type = str)
  parser.add_argument('-D', help = "Required:  SNV table", type = str)
  parser.add_argument('-overlap', help = "Required:  file contain overlapped SNV and CNA", type = str)
  parser.add_argument('-out', help = "Required: Output file", type = str)
  parser.add_argument('-alpha', help = "False positive value", type = float, default= 0.0174)
  parser.add_argument('-beta', help = "False negative value", type = float, default = 0.1256)
  parser.add_argument('-sigma', help = "Standar deviation for snv proportion and cn proportion", type = float, default = 0.001)
  parser.add_argument("-searches", help="Number of searches", type=int, default=100000)
  parser.add_argument("-reveal", help = "Reveal file: where SNV should be placed", type  = str)
  parser.add_argument("-restart", type = int)
  parser.add_argument("-fixCell", help="Whether the cell placement is fixed", type = int, default = 1)
  args = parser.parse_args()

    #start_time = time.time()
  simulate_main(args.tree, args.D, args.overlap, args.reveal, args.out, args.searches, args.alpha, args.beta, args.sigma,  args.restart, args.fixCell)
    #print("--- %s seconds ---" % (time.time() - start_time))
