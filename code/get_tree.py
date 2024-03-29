# tree structure

from copy import deepcopy
import math
from write_tree import *
from TREE import *
import os.path


# this function read the eddgelist
def readTreeStructure(treeFile):
  tree = TREE()
  treeF = open(treeFile, "r")
  lines = treeF.readlines()
  nodes = []
  parent2child = {} # dictionary {parent: [child1, child2]}
  for line in lines:
    temp = line.split()
    nodeID = int(temp[0])
    parentID = int(temp[1])
    if parentID not in parent2child:
      parent2child[parentID] = [nodeID]
    else:
      parent2child[parentID].append(nodeID)
    edgelen = float(temp[3])
    perc = float(temp[4])
    ifLeaf = int(temp[5])
    if ifLeaf:
      tree.leavesID.append(nodeID)
    nodeName = nodeID
    if len(temp) == 7:
      nodeName = temp[6]
    node = MyNode1(id = nodeID, name = nodeName, parent=parentID, parentID=parentID, perc = perc, if_leaf=ifLeaf)
    node.edge_length = edgelen
    nodes.append(node)
  for node in nodes:
    key = node.id
    if key in parent2child:
      node.children = parent2child[key]
  nodes.sort(key = sortNodesKey)
  tree.nodes = deepcopy(nodes)
  return tree

def readD(Dfile):
  f = open(Dfile, "r")
  lines = f.readlines()
  D = []
  for line in lines:
    temp = line.split()
    d = []
    for i in temp:
      val = 0
      if i == "3":
        val = -1
      else:
        val = int(i)
      d.append(val)
    D.append(d)
  return D



def readNewSimulatedTree(treeFile, DFile, overlapped, revealFile, initAlpha, initBeta, initSigma, SNVFile = None, SNVcell =None):
  tree = readTreeStructure(treeFile)
  # read D
  D = readD(DFile)
  tree.D = deepcopy(D)
  SNVnames = [] # list of SNV name
  if SNVFile:
    if os.path.exists(SNVFile):
      f = open(SNVFile)
      lines = f.readlines()
      for line in lines:
        temp = line.rstrip().split()
        name = temp[1].split("_")[0]
        SNVnames.append(name)
      f.close()
  snvs = []
  for m in range(len(tree.D[0])):
    if len(SNVnames) > 0:
      snv = MySNV1(id = m, name = SNVnames[m])
      snvs.append(snv)
    else:
      snv = MySNV1(id = m)
      snvs.append(snv)
  # get SNV cell names
  SNVcellNames = []
  if SNVcell:
    if os.path.exists(SNVcell):
      f = open(SNVcell)
      lines = f.readlines()
      for line in lines:
        temp = line.rstrip().split()[1]
        SNVcellNames.append(temp)
      f.close()
  # get overlapping SNV and CNA
  if os.path.exists(overlapped):
    f = open(overlapped)
    lines = f.readlines()
    for line in lines:
      temp = line.split()
      node = int(temp[0])
      snv = int(temp[1])
      snvs[snv].loss_node.append(node)
    f.close()

  tree.snvs = deepcopy(snvs)
  tree.cellNames = SNVcellNames
  # for d in D:
  #   print(d)
  tree.cells = [-1] * len(tree.D)
  tree.cells_pos = [[] for _ in range(len(tree.D))]
  #read the reveal of cell placement
  if os.path.exists(revealFile):
    f = open(revealFile, "r")
    lines = f.readlines()
    #print(lines)
    for line in lines:
      #print(line)
      temp = line.split()
      if len(temp) < 2:
        continue
      node = int(temp[0])
      if tree.nodes[node].if_leaf:
        # if leaf
        cells = temp[1].split(";")
        if temp[1] != "NA":
          for c in cells:
            c = int(c)
            if node not in tree.cells_pos[c]:
              tree.cells_pos[c].append(node) # if multiple leaves
            else:
              # if previous internal edge contain this leaf
              # set to this leaf
              tree.cells_pos[c] = [node] 
      else:
        leaves = []
        getLeaves(tree, node, leaves)
        cells = temp[1].split(";")
        for c in cells:
          c = int(c)
          # if not set
          if len(tree.cells_pos[c]) == 0:
            tree.cells_pos[c] = deepcopy(leaves)
          # if read a new node with few leaves for this cell
          elif len(tree.cells_pos[c]) > len(leaves) and len(leaves) > 0:
            tree.cells_pos[c].clear()
            tree.cells_pos[c] = deepcopy(leaves)
    f.close()

  cellPosAll = []
  for i in range(len(tree.cells_pos)):
    if len(tree.cells_pos[i]) > 0:
      for p in tree.cells_pos[i]:
        if p not in cellPosAll:
          cellPosAll.append(p)
  tree.cellPosAll = cellPosAll

  n = len(tree.D)
  m = len(tree.D[0])
  X = 0
  for i in range(n):
    for j in range(m):
      if(tree.D[i][j]) == -1:
        X = X + 1
  initGamma = X / (n * m) 
  if(initGamma == 0):
    initGamma = 0.0001
  tree.theta = [initAlpha, initBeta, initGamma]
  tree.sigma = initSigma
  tree.beta_parameters = [
    [2, math.ceil(2*(1 - initAlpha) / initAlpha)], #alpha
    [2, math.ceil(2*(1 - initBeta) / initBeta)], #beta
    [2, math.ceil(2*(1 - initGamma) / initGamma)], #gamma, mean of initial gamma
    [2, math.ceil(2*(1 - initSigma) / initSigma)]] #sigma, mean of signma
  return tree, initGamma

