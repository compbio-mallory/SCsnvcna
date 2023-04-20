# MCMC movement
from copy import deepcopy
import numpy as np

from scipy import stats
import copy
from write_tree import *

from TREE import *
import random

#   # return True if new theta is accepted
def search_theta(tree, m_h, P_DG, old_P_DMK_list, old_P_CM_list):
  old_theta = deepcopy(tree.theta)
  tree.theta, P_DG = updateTheta(tree, 0, m_h)
  tree.theta, P_DG = updateTheta(tree, 1, m_h)
  tree.theta, P_DG = updateTheta(tree, 2, m_h)
  if old_theta != tree.theta:
    # update P_CM_List
    old_P_CM_list_ = update_PCM_list(tree)

    old_P_CM_sum = sum(old_P_CM_list_)
    tree.p = jointP(tree, P_DG, sum(old_P_DMK_list), old_P_CM_sum)
    return True, P_DG, old_P_DMK_list, old_P_CM_list_
  return False, P_DG, old_P_DMK_list, old_P_CM_list
  
# this function search for new sigma
def search_sigma(tree, m_h, P_DG, old_P_DMK_list, old_P_CM_list):
  old_P_DMK_sum = sum(old_P_DMK_list)
  old_P_CM_sum = sum(old_P_CM_list)

  sigma_new = abs(np.random.normal(tree.sigma, 1, 1)[0])
  while(sigma_new > 1 or sigma_new == 0):
    sigma_new = abs(np.random.normal(tree.sigma, 1, 1)[0])

    a = tree.beta_parameters[3][0]
    b = tree.beta_parameters[3][1]
    new_P_DMK = []
    R_sigma = 1
    for j in range(len(tree.snvs)):
      new_P_DMK.append(np.log(getP_DMK(tree.V_bar[j], tree.V_hat[j], sigma_new)))

    new_P_DMK_sum = sum(new_P_DMK)
  
    prior0 = get_beta(sigma_new, a, b)
    prior1 = get_beta(tree.sigma, a, b)

    R_sigma = (new_P_DMK_sum - old_P_DMK_sum) +  (np.log(prior0) - np.log(prior1))
    if(R_sigma > 0):
      R_sigma = 1
    else:
      R_sigma = np.exp(R_sigma)

    if R_sigma >= 1:
      tree.sigma = sigma_new
      P_DMK_list = new_P_DMK[:]
      old_P_DMK_sum = sum(P_DMK_list)
      tree.p = jointP(tree, P_DG, old_P_DMK_sum, old_P_CM_sum)
      return True
    elif(np.random.uniform(0, 1, 1) < m_h):
      if(np.random.uniform(0, 1, 1) < R_sigma):
        tree.sigma = sigma_new
        P_DMK_list = copy.deepcopy(new_P_DMK)
        old_P_DMK_sum = sum(P_DMK_list)
        tree.p = jointP(tree, P_DG,  old_P_DMK_sum, old_P_CM_sum)
        return True
# this function search for new SNV placement
def search_SNV(tree, m_h, P_DG, old_P_DMK_list, old_P_CM_list, ifFixed):
  # print("---------------------")
  new_nodes = deepcopy(tree.nodes)
  new_snvs = deepcopy(tree.snvs)
  new_cells = deepcopy(tree.cells)
  s = np.random.randint(len(tree.snvs), size = 1)
  s = s[0]
  snv = new_snvs[s]
  old_edge = snv.edge
  #print(old_edge)
  edge_list = list(range(len(tree.nodes)))
  # print("edge list", edge_list, "old_edge", old_edge)
  edge_list.remove(old_edge)
  edge_new = random.sample(edge_list,  1)
  edge_new = edge_new[0]
  while edge_new == old_edge:
    edge_new = random.sample(edge_list,  1)
    edge_new = edge_new[0]
  new_snvs[s].edge = edge_new
  removeSNV(new_nodes, new_nodes[old_edge], snv.id)
  new_snvs[s].loss = -1
  new_V_hat = copy.deepcopy(tree.V_hat)
  new_nodes[edge_new].new_snvs.append(snv.id)
  addSNV(new_nodes, new_nodes[edge_new], snv.id)

      
  # handle snv loss due to copy number loss
  new_V_hat[s] = new_nodes[edge_new].perc
  #ifLoss = False
  if len(new_snvs[s].loss_node) > 0:
    subtree = []
    getSubtree(tree.nodes, tree.nodes[edge_new], subtree)
    loss_nodes = [] # snvs can be lost in multiple places
    for loss_n in new_snvs[s].loss_node:
      if loss_n in subtree:
        loss_nodes.append(loss_n)
    loss_node = -1
    if len(loss_nodes) >= 1: # randomly loss on one edge
      loss_node = random.sample(loss_nodes, 1)[0]
    if loss_node in subtree:
      loss = np.random.uniform(0, 1, 1)
      if loss >= 0.5:
        #ifLoss = True
        #print("mutation loss", snv.id, loss_node, new_nodes[edge_new].perc, new_nodes[loss_node].perc)
        new_V_hat[s] = new_nodes[edge_new].perc - new_nodes[loss_node].perc
        removeSNV(new_nodes, new_nodes[loss_node], snv.id)
        if loss_node == edge_new:
          new_nodes[edge_new].new_snvs.append(snv.id)
        new_snvs[s].loss = loss_node
        if new_nodes[loss_node].loss_snvs: 
          new_nodes[loss_node].loss_snvs.append(s)
        else:
          new_nodes[loss_node].loss_snvs = [s]

  cell_move = np.random.uniform(0,1,1)[0]
  new_P_CM = deepcopy(old_P_CM_list)
  if cell_move <= 1:
    affectedCells, affectedLeaves = SNVAffectedCells(new_nodes, tree.leavesID, old_edge, edge_new)

    new_P_CM = updateCellsPlacement(new_nodes, tree.cells_pos, tree.cellPosAll, new_snvs, tree.leavesID, tree.D, new_cells, tree.theta, affectedCells, affectedLeaves, new_P_CM, ifFixed)
    
    del affectedCells, affectedLeaves

  new_P_CM_sum = sum(new_P_CM)

  new_P_DMK = copy.deepcopy(old_P_DMK_list)
  new_P_DMK[s] = np.log(getP_DMK(tree.V_bar[s], new_V_hat[s], tree.sigma))
  # when p_dmk is very close to zero and hence have same log probability -743
  # 
  if new_P_DMK[s] == old_P_DMK_list[s]:
    # if new v hat is closer to old v bar, make the new probability higher
    if abs(tree.V_bar[s] - new_V_hat[s]) < abs(tree.V_bar[s] - tree.V_hat[s]):
      new_P_DMK[s] += 10
    else: # else smaller 
      new_P_DMK[s] -= 10

  new_P_DMK_sum = sum(new_P_DMK)

  r_m = np.random.uniform(0, 1, 1)
  old_P_DMK_sum = sum(old_P_DMK_list)
  old_P_CM_sum = sum(old_P_CM_list)
  R_m = (new_P_DMK_sum  + new_P_CM_sum)- (old_P_DMK_sum + old_P_CM_sum) 

  accept = False
  if(R_m > 0):
    R_m = 1
  elif(R_m == 0):
    R_m = 1

  if(R_m >= 1):

    accept = True
  elif r_m[0] < m_h:
    R_m = np.exp((new_P_DMK_sum  + new_P_CM_sum)- (old_P_DMK_sum + old_P_CM_sum) )
    # print("mh", R_m)
    if (np.random.uniform(0, 1, 1)[0] < R_m):

      accept = True

  if accept:
    #print("Accepted")
    tree.nodes = copy.deepcopy(new_nodes)

    tree.cells = copy.deepcopy(new_cells)

    tree.snvs = copy.deepcopy(new_snvs)
    new_G = getG(new_nodes, new_snvs, new_cells)
    new_V_bar = getV_bar(new_G)
    tree.G = deepcopy(new_G)

    tree.V_bar = deepcopy(new_V_bar)
    tree.V_hat = deepcopy(new_V_hat)

    P_DG = getP_DG(tree.theta, tree.D, tree.G)
    P_DMK_list = []
    P_DMK_sum = 0
    for j in range(len(tree.D[0])):#for each snv
      p_dm = np.log(getP_DMK(tree.V_bar[j], tree.V_hat[j], tree.sigma))
      P_DMK_list.append(p_dm)
      P_DMK_sum += p_dm

    tree.p = jointP(tree, P_DG, P_DMK_sum, new_P_CM_sum)

    return True, P_DG, P_DMK_list, new_P_CM
  else:
    return False,  P_DG, old_P_DMK_list, old_P_CM_list

# MCMC 
# iter: number of iterations
# tree: TREE()
# m_h: probability of m_h movement
# pi: error rate movement
# p_lamda: sigma movement
def MCMC(iter, tree, m_h, pi, p_lamda, ifFixed, burnin):
  tree, P_DMK_list, P_CM_list, P_DG = initTree(tree)
  print(tree.cellPosAll)
  curr_place = []
  for s in tree.snvs:
    curr_place.append(s.edge)
  print("Initial SNV placement:", curr_place)

  print("Initial cells placement:", tree.cells)

  tree.id = 0
  topTreeList = [deepcopy(tree)]
  #print(P_DMK_list)
  #print(P_CM_list)
  i = 1
  while i <= iter:
    #print("--------",i)
    if(i % 5000 == 0):
      print("search %s" % i)
    r = np.random.uniform(0, 1, 1)[0]
    accepted = False
    if(r < pi):
      #print("search theta")
      accepted, P_DG, P_DMK_list, P_CM_list = search_theta(tree, m_h, P_DG, P_DMK_list, P_CM_list)
      #print(tree.theta)
    elif(r < pi + p_lamda):
      accepted = search_sigma(tree, m_h, P_DG, P_DMK_list, P_CM_list)
    else:
      accepted, P_DG, P_DMK_list, P_CM_list = search_SNV(tree, m_h, P_DG, P_DMK_list, P_CM_list, ifFixed)
    if accepted and i >= burnin:
      tree.id = i
      #f_prob.write(str(i) + "\t" + str(round(tree.p[0], 3)) + "\n")
      if len(topTreeList) == 1 and topTreeList[0].id == 0:
        topTreeList.clear()
        topTreeList.append(deepcopy(tree))
      elif tree.p[0] > topTreeList[0].p[0]:
        topTreeList.clear()
        topTreeList.append(deepcopy(tree))
      elif abs(tree.p[0] - topTreeList[0].p[0]) < 0.001:
        topTreeList.append(deepcopy(tree))

    i += 1
  # return dict of inferred tree
  res = {}
  res["G"] = topTreeList[0].G
  res['cells'] = topTreeList[0].cells
  SNVs = {}
  SNVs_loss = {}
  for s in topTreeList[0].snvs:
    SNVs[s.id] = s.edge
    SNVs_loss[s.id] = s.loss
  res['snv'] = SNVs
  res['snv_loss'] = SNVs_loss
  res['p'] = topTreeList[0].p
  res['theta'] = topTreeList[0].theta
  res['id'] = topTreeList[0].id
  res['sigma'] = topTreeList[0].sigma
  res['treeText'] = genTreeText(topTreeList[0])
  #return topTreeList
  #print("returning")
  #print(res)
  return res
