import os, sys
from numpy import roots
from numpy.core.fromnumeric import sort
from numpy.lib.function_base import append
from numpy.lib.shape_base import split
from data_structure import MySNV1, MyNode1, TREE
import copy
from random import sample
import argparse
import numpy as np
from gen_tree_overlappingCNA import *
from CN import CN
from Gen_Ref_Fa import gen_ref, init_ref, write_ref, read_ref, gen_ref_from_tree
def readSNV(file):
    f = open(file, "r")
    SNVs = []
    i = 0
    for line in f:
        temp = line.split()
        name = temp[0]
        chr = temp[1][3:]
        if(chr == "X"):
            chr = 23
        elif(chr == "Y"):
            chr = 24 
        else:
            chr = int(chr)
        pos = int(temp[2])
        snv = MySNV1(i, name, -1, chr, pos, -1, -1)
        SNVs.append(snv)
        i += 1
    return SNVs

def readSNVTable(file):
    f = open(file, "r")
    D = []
    for line in f:
        temp = line.split()
        temp = temp[1:]
        #print(temp)
        snv = []
        for i in temp:
            snv.append(int(i))
        if(len(snv) > 0):
            D.append(snv)
    return D

def readSegs(file):#read segments information and 
    f = open(file, "r")
    segs = []
    line = f.readline()
    header = line
    while(True):
        line = f.readline()
        if not line:
            break
        segs.append(line.split()[:3])
    return segs


#@file log file name from PAUP
#@t target tree number, default is 1
#@seg_num segments number
def readTree(file, seg_num, t = 1):
    f = open(file, "r")
    line = f.readline()
    flag = True
    targetTree = "Tree " + str(t)
    linkage = [] #[node, parent, edge_length]
    cell2node = {}#{cell_name: node_id}
    cn_profiles = {}
    outgroup = -1
    char2num = {
        'A': 10,
        'B': 11,
        'C': 12,
        'D': 13,
        'E': 14,
        'F': 15}
    allNodes = []
    while(flag and line):
        if(line.startswith(targetTree)):
            print("reading " + targetTree)
            while(not line.startswith("-")):
                line = f.readline()
            line = f.readline()
            while(not line.startswith("-")):
                temp = line.split()# reading linkages
                if(len(temp) == 5): #internal node
                    allNodes.append(temp[0])
                    node = int(temp[0])
                    parent = -1
                    if(temp[1] != "root"):
                        parent = int(temp[1])
                    length = int(temp[2])
                    linkage.append([node, parent, length])
                else: #leaf
                    allNodes.append(temp[0])
                    cell = temp[0]
                    if(cell != "diploid"):
                        node = int(temp[1][1:-1]) #get node id
                    else:#if normal diploid cell
                        node = int(temp[1][1:-2])
                        outgroup = node
                        #line = f.readline()
                        #continue
                    parent = -1
                    if(temp[2] != "root"):#set parent
                        parent = int(temp[2])
                    length = int(temp[3])
                    linkage.append([node, parent, length])
                    cell2node[cell] = node
                line = f.readline()   
            for node in allNodes:
                cn_profiles[node] = []
            #reading states
            print(cn_profiles)
            while(not line.startswith("Tree")):
                line = f.readline()
                temp = line.rstrip().split()
                if(len(temp) > 0):
                    start = temp[0]
                    if(start in allNodes):
                        state = temp[1]
                        for i in range(len(state)):
                            cn = state[i]
                            if cn in char2num:
                                cn = char2num[cn]
                            else:
                                cn = int(cn)
                            cn_profiles[start].append(cn)
            flag = False

        line = f.readline()
    return linkage, cell2node, cn_profiles, outgroup

# construct MyNode1
# return a list of MyNode1
#@linkage tree liknage from read_Tree [[node, parent, branch_len], ...]
#@cell2node dict of cell to node id [{"cell": id}]
#@cn_profile dict of copy number profiles for each node {"node": []}
#@segs segments

def makeNode(linkage, cell2node, cn_profiles, segs, outgroup):
    Nodes = []
    cell_list = list(cell2node.keys())
    node_list = list(cell2node.values())
    print(cell_list)
    print(node_list)
    reroot = False
    for i in range(len(linkage)):
        #print("##################################################### i is " + str(i))
        id = linkage[i][0]
        parent = linkage[i][1]
        name = id
        node = MyNode1(id, name, parent)
        node.parentID = int(parent)
        node.parent = parent
        node.edge_length = linkage[i][2]
        if(id in node_list):#update node name
            pos = node_list.index(id)
            name = cell_list[pos]
            node.name = name
            node.if_leaf = True
        if(node.parentID in node_list):#update parent name
            pos = node_list.index(node.parentID)
            node.parent = cell_list[pos]
        #update children
        for j in range(len(linkage)):
            if id == linkage[j][1]:
                node.children.append(linkage[j][0])     
        cn_summary = [{} for _ in range(24)] #
        for j in range(len(segs)):
            chr = int(segs[j][0])
            seg = segs[j][1] + "." + segs[j][2]
            key = str(name)
            val = cn_profiles[key][j]
            cn_summary[chr - 1][seg] = val
            if(node.parentID == -1 and val != 2): #whether to reroot
                reroot = True
        
        node.cn_summary = cn_summary
        #shift id to start with 0
        if(id == outgroup):# if diploid
            node.id = len(linkage) -1
            node.parentID = -1
            node.parent = -1
            child = int(parent) - 2
            node.children.append(child)#change parent to  child, reroot
            node.if_leaf = False
            Nodes.append(node)
            continue
        if(id < outgroup):
            if(id == node.name):
                node.name = node.name - 1    
            node.id = node.id - 1
        else:
            if(id == node.name):
                node.name = node.name - 2   
            node.id = node.id - 2
        #update parent id
        if(node.parentID > -1):
            if(node.parentID < outgroup):
                if(node.parent == node.parentID):#update internal parent name
                    #print("update both parent id and name")
                    node.parent = node.parent - 1
                    node.parentID = node.parentID - 1
            else:  
                if(node.parent == node.parentID):#update internal parent name
                    node.parent = node.parent - 2
                    node.parentID = node.parentID - 2
        for j in range(len(node.children)):
            if (node.children[j] < outgroup):
                node.children[j] =  node.children[j] - 1
            if (node.children[j] > outgroup):
                node.children[j] =  node.children[j] - 2   
               
        Nodes.append(node)
    Nodes.sort(key = sortNodesKey)
    Nodes[len(Nodes) - 2].parentID = len(Nodes) - 1
    Nodes[len(Nodes) - 2].parent = "diploid"
    Nodes[len(Nodes) - 2].edge_length = Nodes[len(Nodes) - 1].edge_length
    Nodes[len(Nodes) - 2].children.remove(outgroup)
    if not reroot:
        Nodes[len(Nodes) - 2].parentID =  - 1
        Nodes[len(Nodes) - 2].parent = -1
        Nodes[len(Nodes) - 2].edge_length = 0
        Nodes = Nodes[:-1]
    Nodes[len(Nodes) - 1].name = "diploid"
    #get percentage of cells in each node
    paths, num_leaf = getPaths(Nodes)
    for i in range(len(Nodes)):
        count = 0
        for j in range(len(paths)):
            if i in paths[j]: 
                count += 1
        Nodes[i].perc = count / num_leaf
    return Nodes
# this function return list of lists path for each leaf to root
# @nodes, list my MyNode1
# return a list of lists, and number of leaves 
def getPaths(Nodes):
    paths = []
    num_leaf = 0
    for node in Nodes:
        if node.if_leaf:
            num_leaf += 1
            curr_node = node
            path = []
            path.append(curr_node.id)
            while curr_node.parentID != -1:
                path.append(curr_node.parentID)
                curr_node = Nodes[curr_node.parentID]
            paths.append(path)
    return paths, num_leaf

def sortNodesKey(node):
    return node.id

# this function check if a snv is overlapped with any cna segments
# @SNVs, list of MySNV1 
# @segs, list of segment from readtree function
# return list of MySNV1 with overlapped segments added
def checkOverlapped(SNVs, segs):
    for i in range(len(SNVs)):
        chrom = SNVs[i].chr
        pos = SNVs[i].pos
        for seg in segs:
            if (chrom == int(seg[0]) and pos >= int(seg[1]) and pos <= int(seg[2])):
                SNVs[i].overlapped = [chrom, int(seg[1]), int(seg[2])]
    return SNVs


def genTreeText(tree):
    root = -1
    nodes = tree.nodes
    numNodes = len(nodes)
    if(nodes[numNodes - 1].parentID == -1):
        root = numNodes - 1
    else:
        for node in nodes:
            if(node.parentID == -1):
                root = node.id
                break
    rootName = nodes[root].name
    rootNHX = ")[&&NHX:N=" + str(rootName) + "];"
    numChild = len(nodes[root].children)
    for c in nodes[root].children:
        child = nodes[c]
        rootNHX = genTreeTextChild(child, nodes) + rootNHX
        numChild -= 1
        if(numChild > 0):
            rootNHX = "," + rootNHX
    rootNHX = "(" + rootNHX +"\n"
    return(rootNHX)

def genTreeTextChild(node, nodes):
    snvs = '+'
    for snv in node.new_snv:
        snvs = snvs + str(snv) + "|"
    snvsLoss = "-"
    for snv in node.loss_snv:
        snvsLoss = snvsLoss + str(snv) + "|"
    cells = '|'
    for cell in node.cells:
        cells = cells + str(cell) + "|"
    currText = ''
    if(node.if_leaf):
        currText = cells + ":"
    e_len = node.edge_length / 5
    if node.edge_length == 0:
        e_len = 0
    currText =  currText  + str(e_len) + "[&&NHX:S="
    if(len(snvs) > 1):
        currText = currText + snvs
    if(len(snvsLoss) > 1):
        currText = currText + ":L=" + snvsLoss
    else:
        currText = currText + ":L="
    cna = "|"
    counter = 0
    for i in range(len(node.cn_summary)):
        for key in node.cn_summary[i]:
            curr_cn = node.cn_summary[i][key] 
            prev = node.parentID 
            prev_cn = 2
            if(key in nodes[prev].cn_summary[i]):
                prev_cn = nodes[prev].cn_summary[i][key]
            if(curr_cn != prev_cn):
                if(curr_cn - prev_cn > 0):
                    cna = cna + "+" + str(counter) + "|"
                if(curr_cn - prev_cn < 0):
                    cna = cna + "-" + str(counter) + "|"
            counter += 1
    currText = currText + ":M=" + cna
    currText = currText +  ":N=" + str(node.name) + "]"
    if(len(node.children) > 0):
        currText = "):" + currText
        numChild = len(node.children)
        for c in node.children:
            child = nodes[c]
            text = genTreeTextChild(child, nodes)
            currText = text + currText
            numChild -= 1
            if(numChild > 0):
                currText = "," + currText
        currText = "(" + currText
    return currText

# this function read the simulated CNA tree 
# @file, path to the file
# return TREE struct
def readSimulatedTree(file):
    tree = np.load(file, allow_pickle=True)
    nodes = []
    segs = []
    oldRootID = 0
    for node in tree:
        if(node.parentID == -1):
            oldRootID = node.id
        newNode = MyNode1(node.id, node.id, parent= node.parentID)
        newNode.parentID = node.parentID

        for key in node.cn_summary:
            if(len(node.cn_summary[key])):
                temp_dict = {}
                for seg in node.cn_summary[key]:
                    temp_dict[seg] = node.cn_summary[key][seg]
                    segSplit = seg.split(".")
                    start = int(segSplit[0])
                    end = int(segSplit[1])
                    chrom = key
                    if([chrom, start, end] not in segs):
                        segs.append([chrom, start, end])
                newNode.cn_summary[key]= temp_dict
        #print(newNode.cn_summary)
        newNode.if_leaf = node.if_leaf
        nodes.append(newNode)
    nodes.sort(key = sortNodesKey)
    #create root
    id = len(nodes) 
    root = MyNode1(id, "root", parent = -1)
    root.parentID = -1
    nodes[oldRootID].parentID = root.id
    nodes.append(root)
    paths, num_leaf = getPaths(nodes)
    #print(paths)
    # print(num_leaf)
    for i in range(len(nodes)):
        count = 0
        for j in range(len(paths)):
            if i in paths[j]: 
                count += 1
        nodes[i].perc = count / num_leaf
    #get children
    for node in nodes:
        for c in nodes:
            if(c.parentID == node.id):
                node.children.append(c.id)
    return nodes, segs
# this function read the simulated SNV table
# return SNV profile, observed data D
def readSimulatedSNV(file, Nodes):
    leaves = []
    snvs = []
    for node in Nodes:
        if node.if_leaf:
            leaves.append(node.id)
    f = open(file, "r")
    line = f.readline()
    d = {}
    for line in f:
        temp = line.rstrip().split()
        key = temp[0] + "." + temp[1] + "." + temp[2] + "." + temp[3] + "." + temp[4]
        cell = int(temp[5])
        if key in d:
            d[key].append(cell)
        else:
            d[key] = [cell]
    counter = 0
    # give new cell id given the leave id
    leaf2cell = {}
    for leaf in leaves: 
        leaf2cell[leaf] = counter
        counter += 1
    D = [[0] * len(d) for _ in range(len(leaves))] 
    #print(D)
    counter = 0
    for key in d:
        keys = key.split(".")
        chrom = int(keys[0])
        pos = int(keys[1])
        newSNV = MySNV1(id = counter, name = -1, ale = int(keys[4]), ori_nuc = keys[2], new_nuc = keys[3], chr = chrom, pos = pos)
        snvs.append(newSNV)
        for leaf in leaves:
            cell = leaf2cell[leaf]
            if(leaf in d[key]):
                D[cell][counter] = 1
        counter += 1

    return snvs, D

# this function randomly create missing entry 
# @D, observed data
# @alpha, FP rate
# @beta, FN rate
# @gamma, missing rate
# return imputed D
def createError(D, alpha, beta, gamma):
    D_cp = copy.deepcopy(D)
    for i in range(len(D)):
        for j in range(len(D[0])):
            p = np.random.uniform(0, 1, 1)[0]
            if(D_cp[i][j] == 0):
                #print("p is " + str(p) + " alpha is " + str(alpha))
                if(p < alpha):
                    #print(str(i) + " " + str(j) + " fp")
                    D_cp[i][j] = 1
                    #print(D_cp[i][j])
                elif(p < alpha + gamma):
                    D_cp[i][j] = -1
            else:
                #print("p is " + str(p) + " beta is " + str(beta))
                if(p < beta):
                    D_cp[i][j] = 0
                elif(p < beta + gamma):
                    D_cp[i][j] = -1
    
    return D_cp

#this function get all possible snv placements for every snv
# use for evluation 
# @nodes MyNode1. list of ground truth node
# @D_bar matrix of ground truth D_bar
# return list of possible edges for each snv to be placed on
def getSNVPlacement(nodes, D_bar):
    res = []
    for j in range(len(D_bar)):
        temp = []
        for node in nodes:
            if(D_bar[j] == node.perc):
                temp.append(node.id)
        res.append(temp)
    return res

# this function is used to read the .out file from SNV_CNA.py
# return true G and imputed D from a run

def readGD(file):
    f = open(file, "r")
    temp = f.readline()
    n = 0
    m = 0
    temp = f.readline()
    temp = temp.rstrip().split()
    n = int(temp[0])
    m = int(temp[1])
    G = []
    D = []
    for i in range(n):
        temp = f.readline()
        temp = temp.rstrip().split()
        g = []
        for j in range(len(temp)):
            g.append(int(temp[j]))
        G.append(g)
    temp = f.readline()
    for i in range(n):
        temp = f.readline()
        temp = temp.rstrip().split()
        d = []
        for j in range(len(temp)):
            d.append(int(temp[j]))
        D.append(d)
    return G, D

# this function create erros on CNA tree
# @tree TREE() object
# @sd standard deviation to change the percentage of cells on a node
# @return return the modified tree 
def createErrorTree(tree, sd):
    #which node to change first
    nodesLen = len(tree.nodes) 
    print(nodesLen)
    edge = np.random.randint(len(tree.nodes) - 1, size = 1)[0]
    while(edge == 0 or edge == len(tree.nodes) - 1):
        edge = np.random.randint(len(tree.nodes) - 1, size = 1)[0]
    parent = tree.nodes[edge].parentID #get parent ID
    oldSibilingTotal = tree.nodes[parent].perc - tree.nodes[edge].perc
    # update perc for current node and its children
    percToChange = np.random.normal(0, sd, 1)[0]
    # edge = 3
    # parent = 1
    # percToChange = 0.03
    newPerc = tree.nodes[edge].perc + percToChange
    while(newPerc <= 0 or newPerc >= 1 or newPerc > tree.nodes[parent].perc ):
        percToChange = np.random.normal(0, sd, 1)[0]
        newPerc = tree.nodes[edge].perc + percToChange
    oldPerc = tree.nodes[edge].perc #old percentage for node being first changed
    tree.nodes[edge].perc = newPerc
    print("node id " + str(tree.nodes[edge].id) + " is changed by  "  + str(percToChange))
    for i in range(len(tree.nodes[edge].children)): #update children's percentage
        updatePerc(tree, edge)
    # update perc for sibilings
    parentOldPerc = tree.nodes[parent].perc
    tree.nodes[parent].perc = 0
    increase = False
    if(percToChange < 0):
        increase = True
    percToChange = abs(percToChange)
    print("current node is " + str(edge) + " parent is " + str(parent) + " update sibilings")
    updateSibliing(tree, increase, edge, parent, percToChange, oldPerc)
    currentNode = parent
    newParent = tree.nodes[currentNode].parentID
    print("update parents")
    while(newParent != -1):
        if(newParent == nodesLen - 1):
            break
        flag = False
        if(tree.nodes[currentNode].perc < parentOldPerc):
            flag = True
        newParentOldPerc = tree.nodes[newParent].perc
        updateSibliing(tree, flag, currentNode, newParent, percToChange, oldPerc)
        parentOldPerc = newParentOldPerc
        currentNode = newParent
        newParent = tree.nodes[currentNode].parentID


# this function update percentage is one level above
def updateSibliing(tree, increase, currentNode, newParent, percToChange, oldPerc):
    tree.nodes[newParent].perc = 0

    for i in range(len(tree.nodes[newParent].children)):
        childID = tree.nodes[newParent].children[i]
        print("child ID is " + str(childID))
        if childID == currentNode:
            tree.nodes[newParent].perc += tree.nodes[childID].perc
            continue
        percToChange1 = percToChange/(1 - oldPerc)*tree.nodes[childID].perc
        print("perc to chagne is " + str(percToChange1))
        if(increase):     
            tree.nodes[childID].perc += percToChange1
        else:
            tree.nodes[childID].perc -= percToChange1
        print(str(childID) + " new perc is " + str(tree.nodes[childID].perc))
        tree.nodes[newParent].perc += tree.nodes[childID].perc
        updatePerc(tree, childID)#update  children
    

# this is a helper function to update the percentage of cells for each child node
# @tree TREE() object
# @return new tree 

def updatePerc(tree, nodeID):
    totalPerc = tree.nodes[nodeID].perc
    totalOldPerc = 0
    for i in range(len(tree.nodes[nodeID].children)):
        childID = tree.nodes[nodeID].children[i]
        totalOldPerc += tree.nodes[childID].perc
    for i in range(len(tree.nodes[nodeID].children)):
        childID = tree.nodes[nodeID].children[i]
        tree.nodes[childID].perc = totalPerc * tree.nodes[childID].perc / totalOldPerc
        for j in range(len(tree.nodes[childID].children)):
            updatePerc(tree, childID)

