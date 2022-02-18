from os import altsep, name
from typing import Counter
import numpy as np
import sys
import re
from scipy import stats
import copy
import argparse
from numpy.core.fromnumeric import size
from numpy.lib.function_base import append, place
from numpy.lib.shape_base import apply_along_axis
from numpy.lib.type_check import _real_dispatcher
from numpy.lib.ufunclike import isposinf
from numpy.random import beta, gamma, sample
import math, statistics
from read_files import *
from data_structure import *

### function to get P(D|G)
def getP_DG(theta, D, G):
    alpha = theta[0]
    beta = theta[1]
    gamma = theta[2]
    if(D == -1 and G == 0):
        return gamma
    elif(D == -1 and G == 1):
        return gamma
    elif(D == 0 and G == 0):
        return 1 - beta - gamma
    elif(D == 0 and G == 1):
        return alpha
    elif(D == 1 and G == 0):
        return beta
    else:
        return 1 - alpha - gamma 


#function to calculate the D_bar
def getD_bar(G):
    n = len(G)
    m = len(G[0])
    D_bar = [0] * m
    for j in range(m):
        total = 0
        for i in range(n):
            total = total + G[i][j]
        D_bar[j] = total / n
    return D_bar
# D_bar = getP_DG(G)

# this function add SNVs to subtree given the tree and current node
# node: MyNode1; current node
# snv: snv id; snv to be added 
def addSNV(tree, node, snv):
    for i in node.children:
        addSNV(tree, tree.nodes[i], snv)
    node.snvs.append(snv)
# this function remove snv from substree
def removeSNV(tree, node, snv_id):
    #print(node.id)
    if node.id != -1:
        for i in node.children:
            removeSNV(tree, tree.nodes[i], snv_id)
        if(snv_id in node.snvs):
            node.snvs.remove(snv_id)
           # print("Removed snv")
        if(snv_id in node.new_snv):
            node.new_snv.remove(snv_id)
           # print("Removed new snv")
        if(snv_id in node.loss_snv):
            node.loss_snv.remove(snv_id)


# check if a snv is on the path to the root given a leaf
# return 1 if true
def pathToRoot(node, tree, snv):
    curr_node = node 
    flag = 0
    while curr_node.parentID != -1:
        #print("DEBUG snv edge " + str(snv.edge))
        if(snv.edge == curr_node.id):
            flag = 1
        if(snv.id in curr_node.loss_snv):
            flag = 0
            return flag
        curr_node = tree.nodes[curr_node.parentID]
    
    return flag

# P_CM = {
#     "-1.1": gamma,
#     "-1.0": gamma,
#     "1.1":1 - beta - gamma,
#     "1.0": alpha,
#     "0.1": beta,
#     "0.0": 1- alpha - gamma
# }
### function to get P(C|M)
def getP_CM(tree, leaf, snv, theta, D):
    alpha = theta[0]
    beta = theta[1]
    gamma = theta[2]
    isPath = pathToRoot(leaf, tree, snv)
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
#function to sort leaf probability
def leafSortKey(pair):
    return pair[1]

# this function return similarity between cells
# @cell1, [], groud truth data for cell 1
# @cell2, [], ground truth data for cell 2
def getDiff(cell1, cell2):
    res = 0
    for j in range(len(cell1)):
        if(cell1[j] != cell2[j]):
            res += 1
    return res

def placeCellsRandom(tree):
    G = tree.G
    cells = [-1] * len(G)
    leaves = []
    for node in tree.nodes:
        if node.if_leaf:
            leaves.append(node.id)
    for i in range(len(G)):
        leaf = sample(leaves, 1)
        leaves.remove(leaf[0])
        cells[i] = leaf[0]
    tree.cells = copy.deepcopy(cells)



# function to place cells on leaves
def placeCells(tree):
    D = tree.D
    G = tree.G
    theta = tree.theta
    snvs = tree.snvs
    cells = [-1] * len(G)

    for i in range(len(G)):
        #for each cell

        res_p = -1 #result product of probability
        res_l = -1 #leaf id
        pos = [] #leaf, p
        for leaf in tree.nodes: #for each leaf
            if not leaf.if_leaf:
                continue
            prod = 1

            for j in range(len(G[0])): #for each mutation
                p_cym = getP_CM(tree, leaf, snvs[j], theta, D[i][j])
                prod = prod * p_cym

            temp = [leaf.id, prod]
            pos.append(temp)


        pos.sort(key = leafSortKey, reverse = True)
        #print(pos)
        topLeaves = []
        topLeaves.append(pos[0])
        for l in range(1, len(pos)):
            if(pos[l][1] == topLeaves[0][1]):
                topLeaves.append(pos[l])
        if(len(topLeaves) > 1):
            for l in range(len(topLeaves)):
                leafID = topLeaves[l][0]
                if(len(tree.nodes[leafID].cells) > 0):
                    diff = 0
                    for c in tree.nodes[leafID].cells:
                        diff += getDiff(G[i], G[c])
                    if(diff > 0):
                        topLeaves[l][1] = topLeaves[l][1] / diff
        topLeaves.sort(key = leafSortKey, reverse = True)

        res_l = topLeaves[0][0]
        tree.nodes[res_l].cells.append(i)
        cells[i] = res_l
    tree.cells = copy.deepcopy(cells)

## Initial placement of cells
# This functio remove all cells from the leaf
def removeCells(tree):
    for node in tree.nodes:
        if node.if_leaf:
            node.cells.clear()  

#this function check if a snv is overlapped with copy number loss
#tree: list of MyNode1
#node: MyNode1
#snv: MySNV1
#return a list of children nodes that support snv loss due to cn loss
#and a list of parent's copy number and child's copy number
def getOverlapped(tree, node, snv):
    res = []
    if len(snv.overlapped) > 0:
        seg = str(snv.overlapped[1]) + "." + str(snv.overlapped[2])
        #print(snv.chr)
        curr_cn = 2
        if(seg in node.cn_summary[snv.chr - 1]):
            curr_cn = node.cn_summary[snv.chr - 1][seg]
        for c in node.children:
            if(c < 0): print("DEBUG" + str(c))
            child = tree.nodes[c]
            c_cn = 2
            if(seg in child.cn_summary[snv.chr - 1]):
                c_cn = child.cn_summary[snv.chr - 1][seg]
            if c_cn < curr_cn:
                eta = (curr_cn - c_cn) / curr_cn
                temp = [c, eta]
                res.append(temp)
    return res

def getSNVLossNode(tree, node, snv):
    if(len(snv.overlapped) > 0):
        #print("### get snv loss node")
        overlapped = getOverlapped(tree, node, snv)
        #print(overlapped)
        for i in range(len(overlapped)):
            r_eta = np.random.uniform(0, 1, 1)
            nodeID = overlapped[i][0]
            eta = overlapped[i][1]
            if r_eta[0] < eta:
                removeSNV(tree, tree.nodes[nodeID], snv)
                tree.nodes[nodeID].loss_snv.append(snv.id)
                return nodeID
        #print(node.children)
        for child in node.children:
            res = getSNVLossNode(tree, tree.nodes[child], snv)
            if res > -1:
                return res
    return -1
    

#function to calculate P(D_bar｜ M_j^k, sigma)
def getP_DMK(D_bar, D_hat, sigma):
    d_jk = abs(D_hat - D_bar)
    p_d_mk = 1
    if d_jk == 0:
        p_d_mk = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp((-0.5) *(d_jk/sigma) ** 2)
    elif(sigma < 0.001):
        p_d_mk =  (np.exp(0.5) / (2 * np.pi)** (0.5)) 
        p_d_mk = p_d_mk * np.exp(-np.exp(2 * np.log(d_jk) - np.log(2) - 2 * np.log(sigma)))
    else:
        p_d_mk = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp((-0.5) *(d_jk/sigma) ** 2)
    if(p_d_mk < 10**-300):
        p_d_mk = 10**-300
    return p_d_mk

#this function calculate beta prior
def get_beta(x, a, b):
    x1 = x + 0.00001
    x0 = x - 0.00001
    p = stats.beta.cdf(x1, a, b) - stats.beta.cdf(x0, a, b)

    if(p < 10**-200):
        p = 10**-200
 
    return p

#this function update theta 
def updateTheta(D, G, theta, flag, beta_parameters, m_h):
    a = beta_parameters[flag][0]
    b = beta_parameters[flag][1]
    value_old = theta[flag]
    value_new = np.random.normal(value_old, 1, 1)
    if(flag == 0):#sample alpha
        while (value_new[0] + theta[2] > 1) or (value_new[0] <= 0):
            value_new = np.random.normal(value_old, 1, 1)
    if flag == 1: #sample  beta
        while (value_new[0] + theta[2] > 1) or (value_new[0] <= 0):
            value_new = np.random.normal(value_old, 1, 1) 
    if flag == 2: #sample gamma
        while(value_new[0] <= 0 or (value_new[0] + theta[0] > 1) or (value_new[0] + theta[1] > 1)):
            value_new = np.random.normal(value_old, 1, 1)
    value_new = value_new[0]
    if(value_new == 0):
        print("DEBUG value_new error " + str(flag))
    theta_new = copy.deepcopy(theta)
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
        for i in range(len(D)):
            for j in range(len(D[0])):
                if( getP_DG(theta_new, D[i][j], G[i][j]) <= 0):
                    print("DEBUG value_new error " + str(theta_new[0]) + " " + str(theta_new[1]) + " " +str(theta_new[2]))
                temp0 += math.log( getP_DG(theta_new, D[i][j], G[i][j]))
                temp1 += math.log( getP_DG(theta, D[i][j], G[i][j])) 
        R_theta = (temp0 - temp1) + (np.log(prior0) - np.log( prior1))
    if(R_theta >= 0):
        R_theta = 1
    else:
        R_theta = np.exp(R_theta)
    if(np.random.uniform(0, 1, 1) < m_h):
        if(R_theta >= 1):
            return theta_new
        elif(np.random.uniform(0, 1, 1) <= R_theta):
            return theta_new
    elif R_theta >= 1:
        return theta_new
    return theta


def placeSNV(tree):
    D_bar = tree.D_bar
    snvs = tree.snvs
    K = len(tree.nodes) - 1  #doesn't place on root
    P_M = 1/K
    sigma = tree.sigma
    D_hat = [0] * len(D_bar)
    #print("PM is " + str(P_M))
    for j in range(len(D_bar)):
        pos = [] #edge, prod, tree, d_hat
        for k in range(K):
            tree1 = copy.deepcopy(tree)
            tree1.nodes[k].new_snv.append(j)
            addSNV(tree1, tree1.nodes[k], j)#add snv
            d_hat = tree1.nodes[k].perc
            loss = getSNVLossNode(tree1, tree1.nodes[k], snvs[j]) #snv loss
            if(loss > -1): d_hat = tree1.nodes[k].perc - tree1.nodes[loss].perc
            p_dmk = getP_DMK(D_bar[j], d_hat, sigma)
            #print("p_dmk is " + str(p_dmk))
            total = 0
            for kk in range(K):
                d_hat1 = tree1.nodes[kk].perc
                if(kk == k): d_hat1 = d_hat
                #print(a)
                total += getP_DMK(D_bar[j], d_hat1, sigma) * P_M
            prod = 0
            if(total != 0.0):
                prod = p_dmk * P_M / (total)
            temp = [k, prod, tree1, d_hat]
            #if(loss > -1):
            #    print("snv loss supported " + str(prod))
            pos.append(temp)
        pos.sort(key = leafSortKey, reverse = True)
        topEdges = []
        topEdges.append(pos[0])
        for e in range(1, len(pos)):
            if(pos[e][1] == topEdges[0][1]):
                topEdges.append(pos[e])
        x = np.random.randint(len(topEdges), size = 1)
        x = x[0]
        tree = copy.deepcopy(topEdges[x][2])
        snvs[j].edge = pos[x][0]
        D_hat[j] = topEdges[x][3]
        if(D_hat[j] != tree.nodes[pos[x][0]].perc):
            print("snv loss supported")
    tree.snvs = copy.deepcopy(snvs)
    tree.D_hat = copy.deepcopy(D_hat)
    return tree

def placeSNVRandom(tree):
    D_bar = tree.D_bar
    K = len(tree.nodes) - 1 #doesn't place on root
    snvs = tree.snvs
    D_hat = tree.D_hat
    D_hat = [0] * len(D_bar)
    #K = len(tree.nodes) - 1
    for j in range(len(D_bar)):
        k = np.random.randint(K - 1, size = 1)
        res_k = k[0]
        addSNV(tree, tree.nodes[res_k], j)
        D_hat[j] = tree.nodes[res_k].perc #initial D_hat
        snvs[j].edge = res_k  
    tree.D_hat = copy.deepcopy(D_hat)
    tree.snvs = copy.deepcopy(snvs)  
    return tree

#function to initialize theta, P_G, P_G0, G, D_bar, D_hat, snv placement, cell placement etc
def init(D, theta, sigma, tree, snvs, beta_parameters):
    #init theta
    X = 0
    num_snv = 0
    n = len(D)
    m = len(D[0])
    for i in range(n):
        for j in range(m):
            if(D[i][j]) == -1:
                X = X + 1
            if(D[i][j] == 1):
                num_snv = num_snv + 1

    #init P_G
    num_non_snv = n * m - X - num_snv
    P_G = (num_snv - num_snv * theta[0] + num_non_snv * theta[1]) / (num_non_snv + num_snv)
    #P_G0 = 1 - P_G
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
    #print(P_G)
    #init G
    G = [[None]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            key0 = ".".join([str(D[i][j]), "0"])
            key1 = ".".join([str(D[i][j]), "1"])
            #P_G = P_G_list[j]
            P_G0 = 1 - P_G
            p_g0 = (getP_DG(theta, D[i][j], 0) * P_G0) / (getP_DG(theta, D[i][j], 0) * P_G0 + 
            getP_DG(theta, D[i][j], 1) * P_G)
            p_g1 = (getP_DG(theta, D[i][j], 1) * P_G) / (getP_DG(theta, D[i][j], 0) * P_G0 + 
            getP_DG(theta, D[i][j], 1) * P_G)
            x = np.random.uniform(0, 1, 1)[0]
            if(x < p_g1):
                G[i][j] = 1
            else:
                G[i][j] = 0
    # init D_bar
    D_bar = getD_bar(G)
    #init random snv placement
    K = len(tree.nodes) - 1
    #print(K)
    P_M = 1/K
    trees = [None] * 11
    for i in range(11):#
        trees[i] = copy.deepcopy(tree)
        trees[i].snvs = copy.deepcopy(snvs)
        trees[i].G = copy.deepcopy(G)
        trees[i].D = copy.deepcopy(D)
        trees[i].sigma = copy.deepcopy(sigma)
        trees[i].theta = copy.deepcopy(theta)
        trees[i].D_bar = copy.deepcopy(D_bar)
    for i in range(len(trees) - 1):
        trees[i]=  placeSNVRandom(trees[i])
        placeCells(trees[i])
        #trees[i].G = getG(trees[i])
        p = jointP(trees[i], beta_parameters)
        trees[i].p = p

    trees[10]  = placeSNV(trees[10])
    placeCells(trees[10])
    #trees[10].G = getG(trees[10])
    p = jointP(trees[10], beta_parameters)
    trees[10].p = p
    print(trees[10].p)
    trees.sort(key=get_Tree_P)

    return  P_G, P_G0, trees[10]
# function to calculate the probability of a model given the placement of cells
def getP(tree, D, G, D_bar, D_hat, theta, sigma, P_M, cells, snvs):
    P = 0
    for j in range(len(D[0])):#for each snv
        p_dm = getP_DMK(D_bar[j], D_hat[j], sigma)
        p_cm = 1
        for i in range(len(D)):
            temp = getP_CM(tree, tree.nodes[cells[i]], snvs[j], theta, D[i][j])
            p_cm *= (temp)

        top = p_dm * P_M * p_cm
        #print("j is " + str(j) + " top is " + str(top) + " pcm is " + str(p_cm) + " p_dm is " + str(p_dm))
        if top < 10**-300:
            #print("flag " + str(10**-300))
            top = 10**-300
        K = len(tree.nodes) - 1
        total = 0
        pdms = []
        pcms = []
        #print("j is " + str(j) + " top is " + str(top) + " pcm is " + str(p_cm) + " p_dm is " + str(p_dm))
        for x in range(K):#for each edge
            D_hat_temp = tree.nodes[x].perc
            if(snvs[j].edge == x):
                total += top
                continue
            p_dmk = getP_DMK(D_bar[j], D_hat_temp, sigma)
            p_cm1 = 1
            temp_snvs = copy.deepcopy(snvs)
            temp_snvs[j].edge = x
            for i in range(len(D)):
                temp = getP_CM(tree, tree.nodes[cells[i]], temp_snvs[j], theta, D[i][j])
                if(temp <= 0):
                    print("DEBUG temp " + str(temp))
                p_cm1 *= (temp)
            #print("edge is " + str(x) + " top is " + str(top) + " pcm is " + str(p_cm1) + " p_dm is " + str(p_dmk))
            
            temp = p_dmk * p_cm1 * P_M 
            #print("temp is " + str(temp))
            if(temp < 10**-300):
                temp = 10**-300
            total = total + temp
            
            pdms.append(p_dmk)
            pcms.append(p_cm1)
        #print("top is " + str(top))
        temp = top / total

        P += np.log(temp)
        #print("P is " + str(P))
    return P

# this functin get G from the cell placement
def getG(tree):
    snvs = tree.snvs
    cells = tree.cells
    G = [[0] * len(snvs) for _ in range(len(cells))]
    for i in range(len(cells)):
        snvOnLeaf = tree.nodes[cells[i]].snvs
        for j in range(len(snvOnLeaf)):
            k = snvOnLeaf[j]
            G[i][k] = 1
    return G 

# return the error rate calcualted based on imputed D and trueG
# return FP, and FN
def getError(D, G):
    fn = 0
    fp = 0
    num1 = 0
    num0 = 0
    for i in range(len(D)):
        for j in range(len(D[0])):
            if(D[i][j] == 1):
                num1 += 1
                if(G[i][j] == 0):
                    fp += 1
            if(D[i][j] == 0):
                num0 += 1
                if(G[i][j] == 1):
                    fn +=1
    FP = fp / num1
    FN = fn / num0
    return FP, FN


# function to write out tree information
def writeTree(f, tree):
    D_bar = tree.D_bar
    D_hat = tree.D_hat
    theta = tree.theta
    sigma = tree.sigma
    G = tree.G
    for row in G:
        f.write(' '.join([str(a) for a in row]) + "\n")
    for item in theta:
        f.write("%s " % item)
    f.write("\n%s\n" % sigma)
    for item in D_bar:
        f.write("%s " % item)
    f.write("\n")
    for item in D_hat:
        f.write("%s " % item)
    f.write("\n")
# This function evaluate the performance of the algorithm
# based on the number of snv placed correctly
# @snvs, list of MySNV1 object 
# @snvPlacement, ground true snv placement. List from getSNVPlacement()
# return a list [total num of snv, num of snv correctly placed, percentage of ]
def evaluation(snvs, snvPlacement):
    count = 0
    for snv in snvs:
        if(snv.edge in snvPlacement[snv.id]):
            count += 1
    total = len(snvs)
    perc = count / total
    return total, count, perc

def jointP(tree, beta_parameters):
    numLeaves = 0
    D = tree.D
    for node in tree.nodes:
        if(node.if_leaf == True):
            numLeaves += 1
    P_C = 1 / numLeaves
    X = 0
    num_snv = 0
    n = len(D)
    m = len(D[0])
    for i in range(n):
        for j in range(m):
            if(D[i][j]) == -1:
                X = X + 1
            if(D[i][j] == 1):
                num_snv = num_snv + 1
  
    num_non_snv = n * m - X - num_snv
    P_G = (num_snv - num_snv * tree.theta[0] + num_non_snv * tree.theta[1]) / (num_non_snv + num_snv)
    
    P_DGtheta = 0
    for i in range(n):
        for j in range(m):
            P_DGtheta += np.log(getP_DG(tree.theta, tree.D[i][j], tree.G[i][j]))
    P_M = 1 / (len(tree.nodes) - 1)

    P_MGCsigma = (getP(tree, tree.D, tree.G, tree.D_bar, tree.D_hat, tree.theta, tree.sigma, P_M, tree.cells, tree.snvs))
    P_alpha = np.log(get_beta(tree.theta[0], beta_parameters[0][0], beta_parameters[0][1]))
    P_beta = np.log(get_beta(tree.theta[1], beta_parameters[1][0], beta_parameters[1][1]))
    P_gamma = np.log(get_beta(tree.theta[2], beta_parameters[2][0], beta_parameters[2][1]))
    P_sigma = np.log(get_beta(tree.sigma, beta_parameters[3][0], beta_parameters[3][1]))
    P = P_DGtheta + P_MGCsigma + P_alpha + P_beta + P_gamma + np.log(P_C) + P_sigma + np.log(P_G)

    return [P, P_DGtheta, P_MGCsigma, P_alpha, P_beta, P_gamma, np.log(P_C), P_sigma, np.log(P_G)]


def placePerfectCells(tree):
    tree.cells = [-1] * len(tree.G)
    usedK = []
    for i in range(len(tree.G)):
        snvs = []
        for j in range(len(tree.G[i])):
            if(tree.G[i][j] == 1):
                snvs.append(j)
        snvs.sort()
        for k in range(len(tree.nodes) - 1):
            if tree.nodes[k].if_leaf:
                snvs1 = tree.nodes[k].snvs
                snvs1.sort()
                if(snvs == snvs1 and k not in usedK):
                    tree.cells[i] = k
                    usedK.append(k)
                    tree.nodes[k].cells.append(i)
                    break

def main(simulate, cnFile, treeFile, snvFile, snvTable, result, searches, initAlpha = 0.0174, initBeta = 0.1256, initSigma = 0.01, treeNum = 1, imputedFP = -1, imputedFN = -1, imputedMiss = -1, errorSD = 0):  
    D = []
    snvs = []
    nodes = []
    segs = []
    g_nodes = [] #ground true tree structure
    true_G = []
    treeText = result + "_treeText"
    f_text = open(treeText, "w")
    f_log = open(result, "w")

    gtree = TREE()
    true_G_bar = []
    snvPlacement = []
    tree = TREE()
    if(simulate == 1):
        g_nodes, segs = readSimulatedTree(treeFile)
        nodes = copy.deepcopy(g_nodes)
        tree.nodes = copy.deepcopy(nodes)
        snvs, true_G = readSimulatedSNV(snvTable, g_nodes) 
        tree.snvs = copy.deepcopy(snvs)
        D = copy.deepcopy(true_G)
      
        tempD = createError(D, imputedFP, imputedFN, imputedMiss)
        tempFP, tempFN = getError(tempD, true_G)
        print("tempFP is " + str(tempFP) + " tempFN is " + str(tempFN) )
        thresholdFP = 0.3
        thresholdFN = 0.3
        if(imputedFN < 0.1):
            thresholdFN = 1
        if(imputedFP < 0.1):
            thresholdFP = 1

        while(abs(tempFP - imputedFP)/imputedFP > thresholdFP or abs(tempFN - imputedFN)/imputedFN > thresholdFN ):
            tempD = createError(D, imputedFP, imputedFN, imputedMiss)
            tempFP, tempFN = getError(tempD, true_G)
            print("tempFP is " + str(tempFP) + " tempFN is " + str(tempFN) )
        D = copy.deepcopy(tempD)

        simulated_tree = np.load(treeFile, allow_pickle=True)
        for node in simulated_tree:
            if(len(node.snvs) > 0):
                for snv in node.snvs:
                    chr = int(snv.chr)
                    pos = int(snv.ref_pos)
                    #print("chr " + str(chr) + " pos " + str(pos))
                    s_id = -1
                    for s in snvs:
                        if(s.chr == chr and s.pos == pos):
                            s_id = s.id
                            parentID = g_nodes[node.id].parentID
                            if(s_id not in g_nodes[parentID].snvs):
                                g_nodes[node.id].new_snv.append(s_id)
                                s.edge = node.id
                            g_nodes[node.id].snvs.append(s_id)
                            break
        gtree.nodes = g_nodes
        gtree.cells = [-1] * len(true_G)
        trueTheta = [initAlpha, initBeta, imputedMiss]
        gtree.theta = trueTheta
        P_M = 1 / (len(g_nodes) - 1)
        D_bar = getD_bar(true_G)
        gtree.D_bar = D_bar
        gtree.D_hat = D_bar
        gtree.sigma = initSigma
        gtree.snvs = copy.deepcopy(snvs)
        K = len(g_nodes) - 1
    else:
        segs = readSegs(cnFile)
        linkage, cell2node, cn_profiles, outgroup = readTree(treeFile, len(segs), treeNum)
        #print(linkage)
        #print(cell2node)
        nodes = makeNode(linkage, cell2node, cn_profiles, segs, outgroup)
        tree.nodes = nodes
        print(genTreeText(tree))
        D = readSNVTable(snvTable)
        tree.D = D
        snvs = readSNV(snvFile)
    
    snvs = checkOverlapped(snvs, segs)
    tree.snvs = copy.deepcopy(snvs)

    X = 0
    n = len(D)
    m = len(D[0])
    for i in range(n):
        for j in range(m):
            if(D[i][j]) == -1:
                X = X + 1
    initGamma = X / (n * m) 
    if(initGamma == 0):
        initGamma = 0.0001
    theta = [initAlpha, initBeta, initGamma]
    print(theta)
    tree.theta = copy.deepcopy(theta)
    tree.sigma = copy.deepcopy(initSigma)
    beta_parameters = [
    [2, math.ceil(2*(1 - initAlpha) / initAlpha)], #alpha, mean of 0.022
    [2, math.ceil(2*(1 - initBeta) / initBeta)], #beta, mean of 0.044
    [2, math.ceil(2*(1 - initGamma) / initGamma)], #gamma, mean of initial gamma
    [2, math.ceil(2*(1 - initSigma) / initSigma)]] #sigma, mean of signma
    for node in tree.nodes:
        print("node id " + str(node.id) + " perc " + str(node.perc))
    if(errorSD > 0):
        createErrorTree(tree, errorSD)
    print("after modification")
    for node in tree.nodes:
        print("node id " + str(node.id) + " perc " + str(node.perc))
    if(simulate == 1):
        gtree.G = copy.deepcopy(true_G)
        gtree.D = copy.deepcopy(D)
        placePerfectCells(gtree)
        print(gtree.cells)
        #recover true theta
        FP_G, FN_G = getError(D, true_G) 
        if(FP_G <= 0.0):
            FP_G = 0.001
        if(FN_G <= 0.0):
            FN_G = 0.001
        gtree.theta = [FP_G, FN_G, initGamma]
        g_text = genTreeText(gtree)
        print(gtree.theta)
        print(g_text)
        print(jointP(gtree, beta_parameters))
        f_text.write("Ground truth tree\n")
        f_text.write(g_text)
    
        f_log.write("Ground truth SNV data\n")
        f_log.write(str(len(true_G)) + " " + str(len(true_G[0])) + "\n")
        for row in true_G:
            f_log.write(' '.join([str(a) for a in row]) + "\n")
        true_G_bar = getD_bar(true_G)
        snvPlacement = getSNVPlacement(nodes, true_G_bar)
        print(snvPlacement)
        
 
    P_G, P_G0, tree = init(D, theta, initSigma, tree, snvs, beta_parameters)

    ## iteration steps
    pi = 0.2#error rate move
    p_lamda = 0.2# sigma move
    p_g = 0 # ground true data move
    m_h = 0.5 #HM moves 
    n = len(D)
    m = len(D[0])
    K = len(tree.nodes) - 1
    print(beta_parameters)

    search = 1
    #print(beta_parameters)
    P = jointP(tree, beta_parameters)
    f_log.write("D\n")
    for row in D:
        f_log.write(' '.join([str(a) for a in row]) + "\n")
    f_log.write("Starting probability is %s\n" % P)
    writeTree(f_log, tree)
    tree.p = 0
    topTrees = []
    for node in tree.nodes:
        if(node.parentID == -1):
            print("root node is " + str(node.id))
    while search <= searches:
        r = np.random.uniform(0, 1, 1)
        r = r[0]
        tree.id = copy.deepcopy(search)
        if(r < pi):
            #print("iteration " + str(search) + " sample new theta\n")
            #update alpha 
            theta= updateTheta(tree.D, tree.G, tree.theta, 0, beta_parameters, m_h)
            tree.theta = copy.deepcopy(theta)
            #update beta
            theta = updateTheta(tree.D, tree.G, tree.theta, 1, beta_parameters, m_h)
            tree.theta = copy.deepcopy(theta)
            #update gamma
            theta = updateTheta(tree.D, tree.G, tree.theta, 2, beta_parameters, m_h)
            tree.theta = copy.deepcopy(theta)

        elif(r < pi + p_lamda):
            #new sigma   

            sigma_new = np.random.normal(tree.sigma, 1, 1)
            while(sigma_new > 1 or sigma_new <= 0):
                sigma_new = np.random.normal(tree.sigma, 1, 1)
            sigma_new = sigma_new[0]
            a = beta_parameters[3][0]
            b = beta_parameters[3][1]
            temp0 = 0
            temp1 = 0
            R_sigma = 1
            for j in range(m):
                k = snvs[j].edge
                temp0 +=  math.log(getP_DMK(tree.D_bar[j], tree.D_hat[j], sigma_new))
                temp1 +=  math.log(getP_DMK(tree.D_bar[j], tree.D_hat[j], tree.sigma))
            prior0 = get_beta(sigma_new, a, b)
            prior1 = get_beta(tree.sigma, a, b)
            R_sigma = (temp0 - temp1) +  (math.log(prior0) - math.log(prior1))
            if(R_sigma > 0):
                R_sigma = 1
            else:
                R_sigma = math.exp(R_sigma)
            if(np.random.uniform(0, 1, 1) < m_h):

                if(R_sigma >= 1):
                    tree.sigma = copy.deepcopy(sigma_new)
                elif(np.random.uniform(0, 1, 1) < R_sigma):
                    tree.sigma = copy.deepcopy(sigma_new)
            elif R_sigma >= 1:
                tree.sigma = copy.deepcopy(sigma_new)
        else:

            tree_new = copy.deepcopy(tree)
            #printTreeinfo(tree)
            s = np.random.randint(m, size = 1)
            s = s[0]
            snv = tree_new.snvs[s]
            old_edge = snv.edge
            edge_new = np.random.randint(K, size = 1)
            edge_new = edge_new[0]
            while edge_new == old_edge:
                edge_new = np.random.randint(K, size = 1)
                edge_new = edge_new[0]
   
            tree_new.snvs[s].edge = edge_new
            removeSNV(tree_new, tree_new.nodes[old_edge], snv.id)
  
            addSNV(tree_new, tree_new.nodes[edge_new], snv.id)
            tree_new.nodes[edge_new].new_snv.append(snv.id)

            tree_new.D_hat[s] = tree_new.nodes[edge_new].perc

            loss = getSNVLossNode(tree_new, tree_new.nodes[edge_new], snv)
            if(loss > -1):
  
                tree_new.D_hat[s] = tree_new.D_hat[s] - tree_new.nodes[loss].perc
            #cells replacement
            removeCells(tree_new)

            placeCells(tree_new)
            tree_new.G = getG(tree_new)

            tree_new.D_bar = getD_bar(tree_new.G)
            #sample r_m
            r_m = np.random.uniform(0, 1, 1)
            R_m = 1
            p_dmk0 = 0
            p_dmk1 = 0
            p_cm0 = 0
            p_cm1 = 0
            for i in range(n):
                p_dmk0 += math.log(getP_DMK(tree_new.D_bar[s], tree_new.D_hat[s], tree_new.sigma))
                p_dmk1 += math.log(getP_DMK(tree.D_bar[s], tree.D_hat[s], tree.sigma))
  
                p_cm0 += math.log(getP_CM(tree_new, tree_new.nodes[tree_new.cells[i]], tree_new.snvs[s], tree_new.theta, tree_new.D[i][s]))
                p_cm1 += math.log(getP_CM(tree, tree.nodes[tree.cells[i]], tree.snvs[s], tree.theta, tree.D[i][s]))

            R_m = (p_dmk0 - p_dmk1) + (p_cm0 - p_cm1) 
            
            if(R_m > 0):
                R_m = 1
            elif(R_m == 0):
                R_m = 0
            else:
                R_m = math.exp(R_m)

            if r_m[0] < m_h:
                #print("mh process: R_m" )
                if R_m >= 1:
                    tree = copy.deepcopy(tree_new)
                elif (np.random.uniform(0, 1, 1) < R_m):
                    tree = copy.deepcopy(tree_new)
            elif(R_m >= 1):
                #print("###### Accept new snv and cell placement3")
                tree = copy.deepcopy(tree_new)

        if (search % 5000 == 0 and search > 0):
            newtree = copy.deepcopy(tree)
            newtree.id = copy.deepcopy(search)
            P = jointP(newtree, beta_parameters)
            newtree.p = P
      
            if(len(topTrees) < 10):
                newtree1 = copy.deepcopy(newtree)
                topTrees.append(newtree1)
                topTrees.sort(key = get_Tree_P)
            elif(newtree.p[0] > topTrees[0].p[0]):
                newtree1 = copy.deepcopy(newtree)
                topTrees[0] = newtree1
                topTrees.sort(key = get_Tree_P)

        search += 1

    for tree in topTrees:
        f_log.write("###### Iteration " + str(tree.id) + "\n")
        f_log.write("Current probability is %s\n" % tree.p[0])
        writeTree(f_log, tree)
    f_log.close()
    for tree in topTrees:
        rootNHX = genTreeText(tree)
        f_text.write(str(tree.id) + "\n")
        f_text.write(rootNHX)
    f_text.close()
    eval = result + "_evaluation"
    f_eval = open(eval, "w")
    if simulate == 1:
        f_eval.write("Tree_ID\tinitialFP\tinitalFN\tinitGamma\tinitialSigma\t" +
        "imputedFP\timputedFN\timputedMiss\testimatedFP\t" +"estimatedFN\testimatedMiss\testimatedSigma\t" +
        "totalSNVs\tcorrectSNVs\tPerc\t" + 
        "FP_D_TrueG\tFN_D_TrueG\tdistFP_D_TrueG\tdistFN_D_TrueG\tdistMiss_G\tProb\t" + 
        "P_DG\tP_MGC\tP_alpha\tP_beta\tP_gamma\tP_C\tP_sigma\tP_G\t" +
        "FP_D_InferredG\tFN_D_InferredG\tdistFP_D_InferredG\tdistFN_D_InferredG\n")
        true_G_bar = getD_bar(true_G)
        print(true_G_bar)
        snvPlacement = getSNVPlacement(nodes, true_G_bar)
        print(snvPlacement)
        trueP = jointP(gtree, beta_parameters)
        f_eval.write(str((-1)) + "\t"  + str(initAlpha) + "\t" + str(initBeta) + "\t" + str(initGamma))
        f_eval.write("\t" + str(initSigma) + "\t" + str(imputedFP) + "\t" + str(imputedFN) + "\t" + str(imputedMiss) + "\t")
        f_eval.write(str(gtree.theta[0]) + "\t" + str(gtree.theta[1]) + "\t" + str(gtree.theta[2]) + "\t")
        f_eval.write(str(round(initSigma, 5)) + "\t" +   str(len(D[0])) + "\t" + str(len(D[0])) + "\t" + str(1) + "\t")
        f_eval.write(str(-1) + "\t" + str(-1) + "\t" + str(-1) + "\t" + str(-1) + "\t-1\t" + str(trueP[0]) + "\t")
        f_eval.write(str(trueP[1]) + "\t" + str(trueP[2]) + "\t" + str(trueP[3]) + "\t")
        f_eval.write(str(trueP[4]) + "\t" + str(trueP[5]) + "\t" + str(trueP[6]) + "\t")
        f_eval.write(str(tree.p[7]) + "\t" + str(tree.p[8]) + "\t" + str(-1) + "\t" + str(-1) + "\t" + str(-1) + "\t" + str(-1) + "\n") 
        for tree in topTrees:
            total, count, perc = evaluation(tree.snvs, snvPlacement)    
            FP_G, FN_G = getError(D, true_G) 
            distFP_D_TrueG = abs(FP_G - tree.theta[0]) 
            distFN_D_TrueG = abs(FN_G - tree.theta[1]) 
            distMiss_G = abs(initGamma - tree.theta[2])
            FP_D_InferredG, FN_D_InferredG = getError(D, tree.G)
            distFP_D_InferredG = abs(FP_D_InferredG - tree.theta[0]) 
            distFN_D_InferredG = abs(FN_D_InferredG - tree.theta[1]) 
            f_eval.write(str((tree.id)) + "\t"  + str(initAlpha) + "\t" + str(initBeta) + "\t" + str(initGamma) + "\t" + str(initSigma))
            f_eval.write("\t" + str(imputedFP) + "\t" + str(imputedFN) + "\t" + str(imputedMiss) + "\t")
            f_eval.write(str(tree.theta[0]) + "\t" + str(tree.theta[1]) + "\t" + str(tree.theta[2]) + "\t")
            f_eval.write(str(tree.sigma) + "\t" +   str(total) + "\t" + str(count) + "\t" + str(perc) + "\t")
            f_eval.write(str(FP_G) + "\t" + str(FN_G) + "\t" + str(distFP_D_TrueG) + "\t" + str(distFN_D_TrueG) + "\t" + str(distMiss_G) + "\t" + str(tree.p[0]) + "\t")
            f_eval.write(str(tree.p[1]) + "\t" + str(tree.p[2]) + "\t" + str(tree.p[3]) + "\t")
            f_eval.write(str(tree.p[4]) + "\t" + str(tree.p[5]) + "\t" + str(tree.p[6]) + "\t")
            f_eval.write(str(tree.p[7]) + "\t" + str(tree.p[8]) + "\t" + str(FP_D_InferredG) + "\t" + str(FN_D_InferredG) + "\t" + str(distFP_D_InferredG) + "\t" + str(distFN_D_InferredG) + "\n") 
    else:
        f_eval.write("Tree_ID\tinitialFP\tinitalFN\tinitGamma\tinitialSigma\t" 
         +"estimatedFP\testimatedFN\testimatedMiss\testimatedSigma\t" + 
        "distFP\tdistFN\tdistMiss\tProb\t" + 
        "P_DG\tP_MGC\tP_alpha\tP_beta\tP_gamma\tP_C\tP_sigma\tP_G\n")
        for tree in topTrees:
            distFP= abs(initAlpha - tree.theta[0]) 
            distFN= abs(initBeta - tree.theta[1]) 
            distMiss= abs(initGamma - tree.theta[2])
            f_eval.write(str((tree.id)) + "\t"  + str(initAlpha) + "\t" + str(initBeta) + "\t" + str(initGamma) + "\t" + str(initSigma) + "\t")
            f_eval.write(str(tree.theta[0]) + "\t" + str(tree.theta[1]) + "\t" + str(tree.theta[2]) + "\t" + str(tree.sigma) + "\t")
            f_eval.write(str(distFP) + "\t" +str(distFN) + "\t" + str(distMiss) + "\t" + str(tree.p[0]) + "\t")
            f_eval.write(str(tree.p[1]) + "\t" + str(tree.p[2]) + "\t" +str(tree.p[3]) + "\t")
            f_eval.write(str(tree.p[4]) + "\t" + str(tree.p[5]) + "\t" +str(tree.p[6]) + "\t")
            f_eval.write(str(tree.p[7]) + "\t" + str(tree.p[8]) + "\n")                  
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cnFile', help = "Required: Merged CN file", type=str)
    parser.add_argument('--treeFile', help = "Required: CN tree log file from PAUP", type = str)
    parser.add_argument('--snvs', help = "Required: SNV file", type = str)
    parser.add_argument('--snvs_table', help = "Required:  SNV table", type = str)
    parser.add_argument('--out', help = "Required: Output file", type = str)
    parser.add_argument('--alpha', help = "False positive value", type = float, default= 0.0174)
    parser.add_argument('--beta', help = "False negative value", type = float, default = 0.1256)
    parser.add_argument('--gamma', help = "Missing rate to create missing simulated data", type = float, default = 0.1)
    parser.add_argument('--sigma', help = "Standar deviation for snv proportion and cn proportion", type = float, default = 0.001)
    parser.add_argument('--treeNum', help = "CNA tree number", type = int, default = 1)
    parser.add_argument('--simulated', help = "Run simulated data", type = int, default = 0)
    parser.add_argument('--imputedFP', help = "FP rate for data imputation", type = float, default = 0.001)
    parser.add_argument('--imputedFN', help = "FN rate for data imputation", type = float, default = 0.001)
    parser.add_argument('--imputedMiss', help = "Missing rate for data imputation", type = float, default = 0.1)
    parser.add_argument("--errorSD", help = "Standar deviation for creating error on CNA tree", type = float, default = 0.0)
    parser.add_argument("--searches", help="Number of searches", type=int, default=1000000)
    args = parser.parse_args()
    print(args)
    if(args.simulated == 1):
        main(args.simulated,  -1, args.treeFile, -1, args.snvs_table, args.out, args.searches, args.alpha, args.beta, args.sigma, -1, args.imputedFP, args.imputedFN, args.imputedMiss, args.errorSD)
    else:
        main(args.simulated, args.cnFile, args.treeFile, args.snvs, args.snvs_table, args.out, args.searches, args.alpha, args.beta,  args.sigma, args.treeNum)
