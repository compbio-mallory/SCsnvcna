# tree structure
from copy import deepcopy
from tkinter.tix import Tree
import numpy as np
from scipy import stats
import copy
from numpy.random import beta, gamma, sample
from write_tree import *
from data_structure import *
import time


class TREE():
  def __init__(self, id = -1, p = 0):
    self.id = id
    self.p = p
    self.D = []
    self.G = []
    self.nodes = []  #list of MYNODE object
    self.V_bar = []
    self.V_hat = []
    self.theta = []
    self.sigma = 0
    self.snvs = [] # list of MYSNV object
    self.cells = [] #where cells are placed
    self.leavesID = [] #list of leave id
    self.beta_parameters = []
    self.cells_pos = None # which nodes should the cell be placed, restricted by the CNA 
    self.cellPosAll = None #list of all reveal leaves position
  
  def get_Tree_P(self):
    return self.p[0]

def sortNodesKey(node):
  return node.id

#function to calculate the V_bar
def getV_bar(G):
  n = len(G)
  m = len(G[0])
  V_bar = [0] * m
  for j in range(m):
    total = 0
    for i in range(n):
      total = total + G[i][j]
    V_bar[j] = total / n
  return V_bar
#this function calculate beta prior
def get_beta(x, a, b):
  x1 = x + 0.00001
  x0 = x - 0.00001
  p = stats.beta.cdf(x1, a, b) - stats.beta.cdf(x0, a, b)
  if(p < 10**-320):
      p = 10**-320
  return p

### function to get P(D|G) for all D and G
def getP_DG(theta, treeD, treeG):
  alpha = theta[0]
  beta = theta[1]
  gamma = theta[2]
  prob = 0
  for i in range(len(treeD)):
    for j in range(len(treeD[i])):
      D = treeD[i][j]
      G = treeG[i][j]
      if(D == -1 and G == 0):
        prob += np.log(gamma)
      elif(D == -1 and G == 1):
        prob += np.log(gamma)
      elif(D == 0 and G == 0):
        prob += np.log(1 - beta - gamma)
      elif(D == 0 and G == 1):
        prob += np.log(beta)
      elif(D == 1 and G == 0):
        prob += np.log(alpha)
      else:
        prob += np.log(1 - alpha - gamma)
  return prob 

### function to get P(D|G) for one point
def getP_DG1(theta, D, G):
  alpha = theta[0]
  beta = theta[1]
  gamma = theta[2]
  prob = 0

  if(D == -1 and G == 0):
    prob += np.log(gamma)
  elif(D == -1 and G == 1):
    prob += np.log(gamma)
  elif(D == 0 and G == 0):
    prob += np.log(1 - beta - gamma)
  elif(D == 0 and G == 1):
    prob += np.log(beta)
  elif(D == 1 and G == 0):
    prob += np.log(alpha)
  else:
    prob += np.log(1 - alpha - gamma)
  return np.exp(prob) 

### function to get P(C|M)
# leaf: Mynode object
# snv: MYSNV1 object
def getP_CM(leaf, snv, theta, D):
  alpha = theta[0]
  beta = theta[1]
  gamma = theta[2]
  #isPath = pathToRoot(leaf, nodes, snv)
  isPath = False
  if snv.id in leaf.snvs:
    isPath = True
  p = 0
  if(D == -1 and isPath == 1):
    p = gamma
  elif(D == -1 and isPath == 0):
    p = gamma
  elif(D == 1 and isPath == 1):
    p = 1 - beta - gamma
  elif(D == 1 and isPath == 0):
    p =  alpha
  elif(D == 0 and isPath == 1):
    p = beta
  else:
    p = 1 - alpha - gamma
  return p

#update log P(C|M) for each cell given the placement of each cell
# due to the change of theta parameters
def update_PCM_list(self):
  P_CM_List = []
  for i in range(len(self.cells)):
    prod = 0
    leafID = self.cells[i]
    for j in range(len(self.G[0])):
      p_cym = getP_CM(self.nodes[leafID], self.snvs[j], self.theta, self.D[i][j])
      prod += np.log(p_cym)
    P_CM_List.append(prod)
  return P_CM_List


def jointP(self, P_DG, P_DMK_sum, P_CM_sum):
  D = self.D
  P_C = 1 / len(self.leavesID)
  num_snv = 0
  n = len(D)
  m = len(D[0])
  X = 0
  for i in range(n):
    for j in range(m):
      if(D[i][j]) == -1:
        X = X + 1
      if(D[i][j] == 1):
        num_snv = num_snv + 1
  num_non_snv = n * m - X - num_snv
  P_G = (num_snv - num_snv * self.theta[0] + num_non_snv * self.theta[1]) / (num_non_snv + num_snv)
  P_DGtheta = P_DG
  P_MG = P_DMK_sum # log p
  P_CM = P_CM_sum
  P_alpha = np.log(get_beta(self.theta[0], self.beta_parameters[0][0], self.beta_parameters[0][1]))
  P_beta = np.log(get_beta(self.theta[1], self.beta_parameters[1][0], self.beta_parameters[1][1]))
  P_gamma = np.log(get_beta(self.theta[2], self.beta_parameters[2][0], self.beta_parameters[2][1]))
  P_sigma = np.log(get_beta(self.sigma, self.beta_parameters[3][0], self.beta_parameters[3][1]))
  P = P_DGtheta + P_MG + P_CM + P_alpha + P_beta + P_gamma  + P_sigma + np.log(P_G)

  return [P, P_DGtheta, P_MG, P_CM, P_alpha, P_beta, P_gamma, P_sigma, np.log(P_G)]


#this function update theta 
def updateTheta(self, flag, m_h):
  a = self.beta_parameters[flag][0]
  b = self.beta_parameters[flag][1]
  value_old = self.theta[flag]
  value_new = abs(np.random.normal(value_old, 0.1, 1)[0])
  if(flag == 0):#sample alpha
    while (value_new + self.theta[2] > 1 or value_new > 0.1) :
      value_new = abs(np.random.normal(value_old, 0.1, 1)[0])
  elif flag == 1: #sample  beta
    while (value_new + self.theta[2] > 1 or value_new > 0.6) :
      value_new = abs(np.random.normal(value_old, 0.1, 1)[0])
  elif flag == 2: #sample gamma
    while((value_new + self.theta[0] > 1) or (value_new + self.theta[1] > 1)):
      value_new = abs(np.random.normal(value_old, 0.1, 1)[0])

  if(value_new == 0):
    print("DEBUG value_new error " + str(flag))
  theta_new = copy.deepcopy(self.theta)
  theta_new[flag] = value_new
  if((theta_new[0] + theta_new[2] > 1) or (theta_new[1] + theta_new[2] > 1)):
    print("DEBUG theta error ")
    print(theta_new)
  temp0 = 0 #new likelihood
  temp1 = 0 #old likelihood
  prior0 = get_beta(value_new, a, b) 
  prior1 = get_beta(value_old, a, b) 
  R_theta = 1
  if(prior0 == 0.0):
    R_theta = 0
  elif(prior1 == 0.0):
    R_theta = 1
  else:
    temp0 = getP_DG(theta_new, self.D, self.G)
    temp1 = getP_DG(self.theta, self.D, self.G)
    R_theta = (temp0 - temp1) +  (np.log(prior0) - np.log( prior1))
    #print("new theta " + str(flag) + " " + str(theta_new[flag]) + " temp0 " + str(temp0) + " temp1 " + str(temp1) + " prior0 " + str(np.log(prior0)) + " prior 1 " + str(np.log(prior1)) )
  if(R_theta >= 0):
    R_theta = 1
  else:
    R_theta = np.exp(R_theta)
  if R_theta >= 1:
    return theta_new, temp0
  elif(np.random.uniform(0, 1, 1) < m_h):
    if(R_theta >= 1):
      return theta_new, temp0
    elif(np.random.uniform(0, 1, 1) <= R_theta):
      return theta_new, temp0
  return self.theta, temp1
  # search for new theta 

def getP_DMK(V_bar, V_hat, sigma):
  d_jk = abs(V_hat - V_bar)
  p_d_mk = stats.norm(0, sigma).pdf(d_jk) 
  if(p_d_mk < 10**-323):
    p_d_mk = 10**-323
  return p_d_mk    

# this function remove SNV from a tree
def removeSNV(nodes, node, snv_id):
  if node.id != -1:
    if(snv_id in node.snvs):
      node.snvs.remove(snv_id)
      
    if(snv_id in node.new_snvs):
      node.new_snvs.remove(snv_id)
      
    if(node.loss_snvs and snv_id in node.loss_snvs):
      node.loss_snvs.remove(snv_id)
    for i in node.children:
      removeSNV(nodes, nodes[i], snv_id)


# this function add an SNV to the tree
def addSNV(nodes, node, snv_id):
  for i in node.children:
    addSNV(nodes, nodes[i], snv_id)
  node.snvs.append(snv_id)

# this function find the cells placement that is affected due to a change of SNV placement
# return a list of cells ID, and a list of leaves ID
def SNVAffectedCells(new_nodes, leavesID, oldEdge, newEdge):
  #get affected leaves ID
  affectedLeaves = []
  affectedCells = []
  subtree = []
  getSubtree(new_nodes, new_nodes[oldEdge], subtree)
  getSubtree(new_nodes, new_nodes[newEdge], subtree)
  #print("DEBUG leave ID", leavesID)
  for leafID in leavesID:
    if leafID in subtree:
      affectedLeaves.append(leafID)
      visited = []
      for cell in new_nodes[leafID].cells:
        if cell not in visited:
          visited.append(cell)
        else:
          print("DEBUG cells on different leave")
          return
        if cell not in affectedCells:
          affectedCells.append(cell)
  return affectedCells, affectedLeaves

#function to sort leaf probability
def leafSortKey(pair):
  return pair[1]

def getSimilarity(c1, c2):
  count = 0
  for i in range(len(c1)):
    if(c1[i] == c1[i]):
      count += 1
  return count
# this funciton place cells on leaves
def placeCells(self):
  D = self.D
  theta = self.theta
  snvs = self.snvs
  cells = [-1] * len(D)
  #print(self.cells_pos)
  for i in range(len(D)):
    #for each cell
    best_leaf = []
    best_leaf_prod = -999999999
    if len(self.cells_pos[i]) > 0:
  
      for leafID in self.cells_pos[i]: #for each leaf
        prod = 0
        for j in range(len(D[0])): #for each mutation
          p_cym = getP_CM(self.nodes[leafID], snvs[j], theta, D[i][j])
          prod += np.log(p_cym)
        if(prod > best_leaf_prod):
          best_leaf.clear()
          best_leaf.append(leafID)
          best_leaf_prod = prod
        elif abs(prod - best_leaf_prod) < 0.000001:
          best_leaf.append(leafID)
          
    else:
      for leafID in self.leavesID: #for each leaf
        if leafID in self.cellPosAll:
          continue
        prod = 0
        for j in range(len(D[0])): #for each mutation
          p_cym = getP_CM(self.nodes[leafID], snvs[j], theta, D[i][j])
          prod += np.log(p_cym)
        if(prod > best_leaf_prod):
          best_leaf.clear()
          best_leaf.append(leafID)
          best_leaf_prod = prod
        elif abs(prod - best_leaf_prod) < 0.000001:
          best_leaf.append(leafID)

    selected = -1
    maxSimilar = -1
    for leaf in best_leaf:
      similarity = 0
      for c in self.nodes[leaf].cells:
        similarity += getSimilarity(D[i], D[c])
      if similarity > maxSimilar:
        maxSimilar = similarity
        selected = leaf
    self.nodes[selected].cells.append(i)
    cells[i] = selected
   
  self.cells = copy.deepcopy(cells)

# this function remove all cells from the leaf
def removeCells(nodes, affectedLeaves):
  for nodeID in affectedLeaves:
    nodes[nodeID].cells.clear()

# this funciton update the cells placement due to the change of SNV placement
def updateCellsPlacement(nodes, cells_pos, cellPosAll, snvs, leavesID, D, cells, theta, affectedCells, affectedLeaves, new_PCM, ifFixed):
  if ifFixed:
    for i in affectedCells:

      prod = 0
      for j in range(len(D[0])): #for each mutation
        p_cym = getP_CM(nodes[cells[i]], snvs[j], theta, D[i][j])
        prod += np.log(p_cym)
      new_PCM[i]=(prod)
    return new_PCM

  # remove the old placement of affected cells
  
  removeCells(nodes, affectedLeaves)
  
  
  for i in affectedCells:

    #for each cell
    best_leaf = []
    best_leaf_prod = -999999999
    #print("old pcm", new_PCM[i], "old placement", cells[i])
    if len(cells_pos[i]) > 0:
      for leafID in cells_pos[i]: #for each leaf
        prod = 0
        for j in range(len(D[0])): #for each mutation
          p_cym = getP_CM(nodes[leafID], snvs[j], theta, D[i][j])
          prod += np.log(p_cym)
        #print("LOH Cell", i, "leaf", leafID, "prod", prod)
        if(prod > best_leaf_prod):
          best_leaf.clear()
          best_leaf.append(leafID)
          best_leaf_prod = prod
        elif abs(prod - best_leaf_prod) < 0.000001:
          best_leaf.append(leafID)
    else:
      for leafID in leavesID: #for each leaf
        if leafID in cellPosAll: # if a cell is not restricted by CNA, then it will not be placed in those leaves
          continue
        prod = 0
        for j in range(len(D[0])): #for each mutation
          p_cym = getP_CM(nodes[leafID], snvs[j], theta, D[i][j])
          prod += np.log(p_cym)
        if(prod > best_leaf_prod):
          best_leaf.clear()
          best_leaf.append(leafID)
          best_leaf_prod = prod
        elif abs(prod - best_leaf_prod) < 0.000001:
          best_leaf.append(leafID)
    selected = -1
    maxSimilar = -1
    for leaf in best_leaf:
      similarity = 0
      for c in nodes[leaf].cells:
        similarity += getSimilarity(D[i], D[c])
      if similarity > maxSimilar:
        maxSimilar = similarity
        selected = leaf
    new_PCM[i]=(best_leaf_prod)
    nodes[selected].cells.append(i)
    cells[i] = selected
  
  return new_PCM

# this function update G due to change of SNV and cell placement
def getG(nodes, snvs, cells):
  G = [[0] * len(snvs) for _ in range(len(cells))]
  for i in range(len(cells)):
    snvOnLeaf = nodes[cells[i]].snvs
    #print("snv on leaf", cells[i], snvOnLeaf)
    for j in range(len(snvOnLeaf)):
      k = snvOnLeaf[j]
      G[i][k] = 1
  return G 

# this funciton place a SNV on the tree
def placeSNV(self):
  V_bar = self.V_bar
  snvs = self.snvs
  K = len(self.nodes) - 1  #doesn't place on root
  sigma = self.sigma
  V_hat = []

  for j in range(len(V_bar)):
    best_prod = -999999
    best_edge = []
    best_perc = []
    for k in range(0, K):
      v_hat = self.nodes[k].perc

      p_dmk = getP_DMK(V_bar[j], v_hat, sigma)
      if p_dmk > best_prod:
        best_prod = p_dmk
        best_edge.clear()
        best_edge.append(k)
        best_perc.clear()
        best_perc.append(v_hat)
      elif p_dmk == best_prod:
        best_edge.append(k)
        best_perc.append(v_hat)
    selected_edge = -1
    selectedP = -1
    diff = 9999
    for e in range(len(best_edge)):
      if abs(best_perc[e] - V_bar[j]) < diff:
        diff =  abs(best_perc[e] - V_bar[j])
        selected_edge = best_edge[e]
        selectedP = best_perc[e]
    snvs[j].edge = selected_edge
    self.nodes[selected_edge].new_snvs.append(j)
    addSNV(self.nodes, self.nodes[selected_edge], j)#add snv

    V_hat.append(selectedP)

  self.V_hat = copy.deepcopy(V_hat)
  


def placeSNVRandom(self):
  V_bar = self.V_bar
  K = len(self.nodes)  #doesn't place on root
  snvs = self.snvs
  V_hat = self.V_hat
  V_hat = [0] * len(V_bar)

  #K = len(self.nodes) - 1
  for j in range(len(V_bar)):
    k = np.random.randint( 1, K ,  1)
    res_k = k[0]
    self.nodes[res_k].new_snvs(j)
    addSNV(self.nodes, self.nodes[res_k], j)
    V_hat[j] = self.nodes[res_k].perc #initial V_hat
    snvs[j].edge = res_k  
  self.V_hat = copy.deepcopy(V_hat)
  self.snvs = copy.deepcopy(snvs)  




def initTree(self):
  D = self.D
  theta = self.theta
  snvs = self.snvs
  nodes = self.nodes
  num_snv = 0
  n = len(D)
  m = len(D[0])
  # print("n", n, "m", m)
  X = 0
  for i in range(n):
    for j in range(m):
      if(D[i][j]) == -1:
        X = X + 1
      if(D[i][j] == 1):
        num_snv = num_snv + 1
  num_non_snv = n *m - X - num_snv
  P_G = (num_snv - num_snv * theta[0] + num_non_snv * theta[1]) / (num_non_snv + num_snv)
  P_G_list = [0] * m
  for j in range(m):
    num_snv = 0
    x = 0
    for i in range(n):
      if(D[i][j] == 1):
        num_snv += 1
      if(D[i][j] == -1):
        x += 1
    non_snv = n - num_snv - x
    pg = (num_snv - num_snv * theta[0] + non_snv * theta[1]) / (non_snv + num_snv)
    P_G_list[j] = pg
  # print("P_G_List", P_G_list)
  #init G
  G = [[None]*m for _ in range(n)]
  for i in range(n):
    for j in range(m):
        #P_G = P_G_list[j]
      P_G0 = 1 - P_G_list[j]
      p_g1 = (getP_DG1(theta, D[i][j], 1) * P_G_list[j]) / (getP_DG1(theta, D[i][j], 0) * P_G0 + 
      getP_DG1(theta, D[i][j], 1) *P_G_list[j])
      x = np.random.uniform(0, 1, 1)[0]
      if(p_g1 >= 0.5):
        G[i][j] = 1
      else:
        G[i][j] = 0
  # init D_bar
  V_bar = getV_bar(G)
  self.V_bar = deepcopy(V_bar)
  self.G = deepcopy(G)

  placeSNV(self)

  placeCells(self)
  P_DMK_list = [] # m muation site
  P_CM_list = [] # n cells 
  for j in range(len(D[0])):#for each snv
    p_dm = np.log(getP_DMK(self.V_bar[j], self.V_hat[j], self.sigma))
    P_DMK_list.append(p_dm)
  
  for i in range(len(D)):
    p_cm = 0
    for j in range(len(D[0])):
      temp = np.log(getP_CM(nodes[self.cells[i]], snvs[j], theta, D[i][j]))
      p_cm += (temp)
    P_CM_list.append(p_cm)
  # print(P_DMK_list)
  # print(P_CM_list)
  P_DG = getP_DG(self.theta, self.D, self.G)
  self.p = jointP(self, P_DG, sum(P_DMK_list), sum(P_CM_list))
  
  return self, P_DMK_list, P_CM_list, P_DG

# This function evaluate the performance of the algorithm
# based on the number of snv placed correctly
# @snvs, list of MySNV1 object 
# @snvPlacement, ground true snv placement. List from getSNVPlacement()
# return a list [total num of snv, num of snv correctly placed, percentage of ]
def evaluation(snvs, snvPlacement):
  count = 0
  for snv in snvs:
    if(snv.edge == snvPlacement[snv.id]):
      count += 1
  total = len(snvs)
  perc = count / total
  return total, count, perc

# this function check the loss of information
# find the SNV cells placements that are restrainted by the 
def LOICNA(tree, cells_cn):
  cells_pos = [[] for _ in range(len(cells_cn))]
  # 
  for i in range(len(cells_cn)):
    cell = cells_cn[i]
    for chr in range(24):
      CNAs1 = cell[chr]
      chr1 = chr
      flag = False
      for key in CNAs1:
        temp = key.split(".")
        start1 = int(temp[0])
        end1 = int(temp[1])
        seg1 = [chr, start1, end1]
        for node in tree.nodes:
          CNAs2 = node.cn_summary[chr]
          for key in CNAs2:
            temp = key.split(".")
            start2 = int(temp[0])
            end2 = int(temp[1])
            seg2 = [chr, start2, end2]
            isSame = ifSameCNA(seg1, seg2, 0.8)
            if isSame:
              leaves = []
              getLeaves(tree, node.id, leaves)
              cells_pos[i] = deepcopy(leaves)
              flag = True
              break
          if flag:
            break
        if flag:
          break
  return cells_pos

# get leaves below a given node
def getLeaves(tree, nodeID, leaves):
  if(not tree.nodes[nodeID].if_leaf):
    for c in tree.nodes[nodeID].children:
      getLeaves(tree, c, leaves)
  if tree.nodes[nodeID].if_leaf:
    leaves.append(nodeID)

# this function check if two CNA segments are overlapped at p%
def ifSameCNA(seg1, seg2, perc = 0.8):
  chr1 = seg1[0]
  chr2 = seg2[0]
  start1 = seg1[1]
  start2 = seg2[1]
  end1 = seg1[2]
  end2 = seg2[2]
  if chr1 != chr2:
    return False
  if end1 <= start2 or end2 <= start1:
    return False

  if start1 <= start2:
    if end1 <= end2:
      overlapped = end1 - start2
      if overlapped / (end1 - start1) >= perc and overlapped / (end2 - start2) >= perc:
        return True
    else:
      overlapped = end2 - start2
      if overlapped / (end1 - start1) >= perc and overlapped / (end2 - start2) >= perc:
        return True
  elif start1 > start2:
    if end2 <= end1:
      overlapped = end2 - start1
      if overlapped / (end1 - start1) >= perc and overlapped / (end2 - start2) >= perc:
        return True
    else:
      overlapped = end1 - start1
      if overlapped / (end1 - start1) >= perc and overlapped / (end2 - start2) >= perc:
        return True
  return False

# this function filter out CNA on tree nodes that are not overlapped with SNV
def filterCNA(tree):
  snvs_pos = []
  for s in tree.snvs:
    if len(s.overlapped) > 0:
      pos = s.overlapped[0]
      if pos not in snvs_pos:
        snvs_pos.append(pos)
  for node in tree.nodes:
    for i in range(24):
      chr = i
      cn_dict = node.cn_summary[i]
      key_to_delete = []
      for key in cn_dict:
        temp = key.split(".")
        start = int(temp[0])
        end = int(temp[1])
        temp1 = [chr, start, end]
        if temp1 not in snvs_pos:
          key_to_delete.append(key)
        else:
          for s in tree.snvs:
            if len(s.overlapped) > 0:
              if temp1 == s.overlapped[0]:
                s.overlapped[1].append(node.id)
      for k in key_to_delete:
        del cn_dict[k]
  # For each overlapped SNVs get the list of involved nodes

# get copy number loss node if place a SNV before the overlapped CNA loss
# @tree: TREE()
# @snv: MYSNV1()
def getSNVLossNode(tree, snv, nodeID):
  flag = False
  resNode = -1
  if(len(snv.overlapped) > 0):
    if nodeID in snv.overlapped[1]:
      parentID = tree.nodes[nodeID].parentID
      flag = ifLoss(tree.nodes[parentID], tree.nodes[nodeID], snv.overlapped[0])
      if flag:
        resNode = nodeID
        return True, resNode
      else:
        for child in tree.nodes[nodeID]:
          flag, resNode = getSNVLossNode(tree, snv, child)
          if flag:
            return flag, resNode
  return flag, resNode

# check if a node has copy number loss
# @parent MyNode1
# #curr_node MyNode1
def ifLoss(parent, curr_node, seg):
  chr = seg[0]
  key = seg[1] + "." + seg[1]
  parentCN = 2
  if key in parent.cn_summary[chr]:
    parentCN = parent.cn_summary[chr][key]
  currCN = curr_node.cn_summary[chr][key]
  if currCN < parentCN:
    return True
  else:
    return False

def getSubtree(nodes, curr_node: MyNode1, subtree):
  if(curr_node.id not in subtree):
    subtree.append(curr_node.id)
  for c in curr_node.children:
    getSubtree(nodes, nodes[c], subtree)
  

def pathToRoot(tree, node):
  path = [node.id]
  curr_node = node.id
  while tree.nodes[curr_node].parentID != -1:
    path.append(tree.nodes[curr_node].parentID)
    curr_node =tree.nodes[curr_node].parentID
  return path