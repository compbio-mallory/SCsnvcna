#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@authors: Xian Fan Mallory
Contacting email: fan@cs.fsu.edu
"""

from ast import Continue
import sys
import argparse
import numpy as np
import random
import copy

class Edge():
    def __init__(self, p, c):
        self.p = p
        self.c = c

class Node():
    def __init__(self, e, p, c):
        # only the edge above this node is listed here
        self.e = e
        self.p = p
        # this node may have multiple children, c is an array
        self.c = c

# add one SNV cell to each internal node, return a new SNVcell_array with the new lines attached to those on the leaves
def add_SNVcells_internal(tree_f, SNVcell_array, cell_n, out_f):

    ret = copy.deepcopy(SNVcell_array)
    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')

    # look for internal nodes, for each internal node, add one cell to it, starting from n, record all the cell attachment (including those on the leaves) to out_f

    p_dict = {}
    while(line != ""):
        line_a = line.split('\t')
        # the 1st (0-based) column is the parent 
        p = line_a[1]
        p_dict[p] = 1
        line = file.readline().rstrip('\n')

    file.close()

    SNVcell_ID = cell_n
    for i in p_dict.keys():
        if i == "-1":
            continue
        ret.append([i, str(SNVcell_ID)])
        SNVcell_ID = SNVcell_ID + 1

    # now add one cell to each leaf node that doesn't have any cell attached to it
    for i in SNVcell_array:
        [leafID, cellIDs] = i
        if cellIDs == "NA":
            # this leaf does not have any cell attached to it, add a cell
            ret.append([leafID, str(SNVcell_ID)])
            ret.remove(i)
            SNVcell_ID = SNVcell_ID + 1 

    print(str(ret))

    print_f(ret, out_f)
    return ret, SNVcell_ID

# initialize a matrix that has r rows and c columns
def init_m(r, c, value):
    m = []
    for i in range(r):
        m.append(copy.deepcopy([value]*c))
    return m



# return a n*m matrix with zero entries
def init(n, m):
    ret = []
    for i in range(n):
        ret.append([0]*m)
        for j in range(m):
            ret[i][j] = 0
    return ret

   
def print_matrix(M, n, matrix_file):
    matrix_f = open(matrix_file, "w")
    for i in range(n):
        str_ = [str(j) for j in M[i]]
        print("\t".join(str_), file = matrix_f) 
    matrix_f.close()

def print_f(array, out_file):
    out_f = open(out_file, 'w')
    for i in range(len(array)):
        print("\t".join(array[i]), file = out_f)
    out_f.close()

# return an array of mutations that are above an edge
def mut_above_edge(edgeID, e_dict, n_dict, mut_array):
    # turn mut_array into a dict. mut_array has edgeID and mut joined by semicolon
    mut_dict = {}
    for i in range(len(mut_array)):
        mut_dict[mut_array[i][0]] = mut_array[i][1]
    p = e_dict[edgeID].p
    ret = []
    if mut_dict[edgeID] != "NA":
        ret = mut_dict[edgeID].split(";")
    while p != "-1":
        edgeID = n_dict[p].e
        if mut_dict[edgeID] != "NA":
            for i in mut_dict[edgeID].split(";"):
                ret.append(i)
        p = e_dict[edgeID].p
    #print(str(ret))
    return ret

# return the array of internal node IDs under this edge
def internal_nodes_under_edge(edgeID, e_dict, n_dict):
    ret_node = []
    node = e_dict[edgeID].c
    # a stack to store all the nodes below for retrieval of the leaves
    c_arr = [node]
    while c_arr != []:
        node = c_arr.pop()
        #print(node)
        if n_dict[node].c != "NA":
            # more children
            for i in n_dict[node].c.split(";"):
                c_arr.append(i)

            # add this internal node to the list
            ret_node.append(node)

    return ret_node
 
# return the array of leaf IDs under this edge
def leaves_under_edge(edgeID, e_dict, n_dict):
    ret_node = []
    node = e_dict[edgeID].c
    # a stack to store all the nodes below for retrieval of the leaves
    c_arr = [node]
    while c_arr != []:
        node = c_arr.pop()
        #print(node)
        if n_dict[node].c != "NA":
            # more children
            for i in n_dict[node].c.split(";"):
                c_arr.append(i)
        else:
            ret_node.append(node)

    return ret_node
        

# given a leaf ID, and an edge dict that has edge ID as the key, .c and .p as the child and parent node ID, return all the edges above this leaf ID in array 
def retrieve_edges(leafID, n_dict):
    p = n_dict[leafID].p
    e = n_dict[leafID].e
    e_array = [e]
    while p != "-1":
        ID = p
        p = n_dict[ID].p
        e = n_dict[ID].e
        e_array.append(e)
    return e_array

# not adding any errors, return the G matrix based on which cell has which mutations
def makeG(n, m, SNVcell_mut): 
    G = init(n, m)
    for cellID in range(n):
        if str(cellID) in SNVcell_mut.keys():
            mut_array_ = SNVcell_mut[str(cellID)]
            for mutID in range(m):
                if str(mutID) in mut_array_:
                    G[cellID][mutID] = 1
    return G

def cell2mutation_arrays(SNVcell_array, node_dict, mut_array_dict):
    # a dict from SNV cell ID to the mutation ID array
    SNVcell_mut = {}

    for i in range(len(SNVcell_array)):
        leafID = SNVcell_array[i][0]
        SNVcellIDs = SNVcell_array[i][1]
        # it's possible no SNV cells is assigned to this leaf
        if SNVcellIDs == "NA":
            continue
        # given the leaf ID, return the edges above this leaf in array
        edges_above_leaf = retrieve_edges(leafID, node_dict)
        # mutation IDs in an array corresponding to this SNV cell
        mutID_array = []
        for e in edges_above_leaf:
            if mut_array_dict[e] != "NA":
                for mutID in mut_array_dict[e].split(";"):
                    mutID_array.append(mutID)
   
        for j in SNVcellIDs.split(";"):
            SNVcell_mut[j] = copy.deepcopy(mutID_array)

    return SNVcell_mut

# copy the first n rows from D_miss to D_miss_, then add the same rate of missing values to the rest of the rows in G_.
def add_missing_SCARLET(D_miss, G_, n, m, missingP):
    D_miss_ = init(len(D_miss), m)

    # copy the first n rows
    for i in range(n):
        for j in range(m):
            D_miss_[i][j] = D_miss[i][j]

    # add missing to the rest of the rows in G_
    G_add = init(len(G_) - n, m)
    # first, copy these rows
    for i in range(len(G_) - n):
        for j in range(m):
            G_add[i][j] = G_[i + n][j]

    # apply add_missing to this matrix 
    G_add_ = add_missing(G_add, len(G_add), m, missingP)

    for i in G_add_:
        D_miss_.append(i)

    return D_miss_

# copy the first n rows from D_miss_FP_FN to D_miss_FP_FN_, for the rest of the rows in D_miss_, apply add_FPFNs
def add_FPFNs_SCARLET(D_miss_FP_FN, D_miss_, n, m, alpha, beta):

    D_miss_FP_FN_ = init(len(D_miss_FP_FN), m)

    # copy the first n rows
    for i in range(n):
        for j in range(m):
            D_miss_FP_FN_[i][j] = D_miss_FP_FN[i][j]

    # add FPs and FNs to the rest of the rows in D_miss_
    D_miss_add = init(len(D_miss_) - n, m)
    # first, copy these rows
    for i in range(len(D_miss_) - n):
        for j in range(m):
            D_miss_add[i][j] = D_miss_[i + n][j]

    total_zeros_ = count_total_value(D_miss_add, len(D_miss_add), m, 0)
    total_ones_ = count_total_value(D_miss_add, len(D_miss_add), m, 1)
    # apply add_FPFNs to this matrix
    D_miss_add_ = add_FPFNs(D_miss_add, len(D_miss_add), m, alpha, beta, total_zeros_, total_ones_)

    for i in D_miss_add_:
        D_miss_FP_FN_.append(i)

    return D_miss_FP_FN_

# given the mutation edge pair (mutation will be lost underneath this edge), the current G matrix, copy the first n rows of G matrix, and add loss to the cells attached to the internal nodes
def addLossMutSCARLET(G_, e_dict, n_dict, SNVcell_array_, n, m, G, mutEdgePair):
    # copy the first rows of G to return G
    retG = init(len(G_), m)
    for i in range(n):
        for j in range(m):
            retG[i][j] = G[i][j]

    # copy the rest of the rows in G_ to return G
    for i in range(len(G_) - n):
        for j in range(m):
            retG[i + n][j] = G_[i + n][j] 

    # according to mutEdgePair, add the losses to those cells attached to an internal node
    # from internal node ID to cell ID for only internal nodes
    # 06022022: also allow the leaf node to appear as some leaf node had "NA" and needs one cell to be attached to it
    dict_node_to_cell = {}
    for i in SNVcell_array_:
        [nodeID, cellID] = i
        #if n_dict[nodeID].c != "NA":
        #    this is an internal node
        #    06022022 can also be a leaf node 
        #    dict_node_to_cell[nodeID] = cellID
        # only deal with those with one cell
        if cellID.isdigit():
            dict_node_to_cell[nodeID] = cellID 

    for i in mutEdgePair:
        [mut, edge] = i
        internalNodes = internal_nodes_under_edge(edge, e_dict, n_dict)
        leafNodes = leaves_under_edge(edge, e_dict, n_dict)
        # mut will be lost on all the cells attached to these internalNodes
        for node in internalNodes + leafNodes:
            if node in dict_node_to_cell.keys():
                c = int(dict_node_to_cell[node])
                if c < n:
                    print("Warning: the cell " + str(c) + " attached to the node " + node + " is within the matrix size " + str(n))
                # mut will be lost on c
                if retG[c][int(mut)] != 1:
                    print("Warning: the cell " + str(c) + " attached to node " + node + " does not have mutation " + mut + ". It might have been changed prior to the SCARLET code since it is the only cell attached to this leaf node. ")
                else:
                    retG[c][int(mut)] = 0

    return retG
          


def addLossMut(G, lossMutP, e_dict, n_dict, mut_array, SNVcell_array, mutLoss_file, mutCNAoverlap_file, m):
    # the edges selected
    edge_array = sorted(e_dict.keys())
    # total edges selected and total mutations lost
    lossN = int(lossMutP * len(edge_array))
    edges = random.sample(range(len(edge_array)), lossN)
    # find all mutations above all edges, use a dict to store them
    mut_cand_dict = {}
    for i in range(len(edges)):
        edge = str(edges[i])
        # find all mutations above all edges, use a dict to store them
        mut_cand = mut_above_edge(edge, e_dict, n_dict, mut_array) 
        for j in mut_cand:
            if j not in mut_cand_dict.keys():
                mut_cand_dict[j] = edge
            else:
                # mutliple edges are below this mutation, record all of them
                mut_cand_dict[j] = mut_cand_dict[j] + ";" + edge

    mut_cand_arr = list(mut_cand_dict.keys())
    # randomly select lossN mutations in the candidate
    if len(mut_cand_arr) < lossN:
        lossN = len(mut_cand_arr)
    muts_array = random.sample(range(len(mut_cand_arr)), lossN)

    # make an array with the mutation IDs that is to be lost 
    muts_array_IDs = []
    for i in range(len(muts_array)):
        muts_array_IDs.append(mut_cand_arr[i])

    mutLoss_f = open(mutLoss_file, "w")
    mutCNAoverlap_f = open(mutCNAoverlap_file, "w")

    # record which mutation is lost after which edge, to be used for SCARLET
    mutEdgePair = []
    # for each selected mutation, randomly select an edge below which the cells will lose this mutation
    for i in range(len(muts_array)):
        mutID = mut_cand_arr[i]
        # edges that have been selected below this mutation
        e_arr = mut_cand_dict[mutID].split(";")
        # randomly select one
        edgeID = e_arr[random.sample(range(len(e_arr)), 1)[0]]
        # find the leaves below this edge
        leafIDs = leaves_under_edge(edgeID, e_dict, n_dict)
        # find all the cell IDs on these leaves
        cellsIDs = []
        for j in range(len(SNVcell_array)):
            if SNVcell_array[j][0] in leafIDs:
                if SNVcell_array[j][1] != "NA":
                    # all those cell IDs on this row of SNVcell_array shall be recorded
                    for k in SNVcell_array[j][1].split(";"):
                        cellsIDs.append(k)

        # change G matrix
        for j in cellsIDs:
            #print(j + "    " + mutID)
            if G[int(j)][int(mutID)] == 0:
                print("Warning: G[" + str(j) + "][" + str(mutID) + "] is zero in loss of mutation operation. ")
            #print("In mutation loss, change cell " + j + " mutation " + mutID + " to 0. ")
            G[int(j)][int(mutID)] = 0

        # record the mutation loss in a file, the first column being the edge ID, the second column being the mutation IDs lost under this edge 
        print("\t".join([edgeID, mutID]), file=mutLoss_f)
        print("\t".join([edgeID, mutID]), file=mutCNAoverlap_f)

        mutEdgePair.append([mutID, edgeID])

    # randomly add more overlaps between SNVs and CNAs (edges on the tree)
    for i in range(m):
        if i not in muts_array_IDs:
            # not the SNVs that have been lost, add overlaps 
            # select edges
            if random.uniform(0, 1) > 0.5:
                # this SNV overlaps with some CNA on the tree
                edgeID = edge_array[random.sample(range(len(edge_array)), 1)[0]]
                print("\t".join([edgeID, str(i)]), file=mutCNAoverlap_f)
                continue
 
    mutLoss_f.close()
    mutCNAoverlap_f.close()
    return G, mutEdgePair

def read_edge_num(tree_f):
    ret = 0
    file = open(tree_f, "r")
    line = file.readline()
    while(line != ""):
        ret = ret + 1
        line = file.readline()
    file.close()
    return ret 

# If no mutation on this edge, then "NA" on the second column
def distribute_mutations(mut_n, tree_f, out_f):
    # to make all mutations add up to mut_n
    edge_num = read_edge_num(tree_f)
    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    mut_array = []
    # mutation ID starts from 0 (0-based)
    prev = -1
 
    # to avoid negative number of mutation, go through the whole file once and if negative, go through each line again
    edge = []
    while(line != ""):
        line_a = line.split('\t')
        # all branch lengths have already been normalized. They add up to 1
        branch_len = float(line_a[3])
        edge.append([line_a[0], branch_len])
        line = file.readline().rstrip('\n')
    file.close()

    tag = True
    mut_nums = init_m(len(edge), 2, 0)
    while tag:
        prev = 0
        for it in range(len(edge)):
            edgeID, branch_len = edge[it]
            # 06132022 Add a protection so that the total will not be above mut_n
            mut_num = round(mut_n * branch_len)

            if prev + mut_num > mut_n or it == len(edge) - 1:
                mut_num = mut_n - prev

            #print("In it " + str(it) + ", mut_num = " + str(mut_num))

            mut_nums[it] = [edgeID, mut_num]

            # 06132022 This warning shall not appear again. 
            if mut_num < 0:
                print("In distributing mutations, mut_num becomes negative. Start again. ")
                break
            prev = prev + mut_num
            if it == len(edge) - 1:
                print("Success in distributing mutations. ")
                tag = False

    prev = -1
    for i in range(len(mut_nums)):
        edgeID, mut_num = mut_nums[i]
        mut_IDs = "NA"
        if mut_num != 0:
            mut_IDs = str(prev + 1)
            for j in range(prev + 2, prev + mut_num + 1):
                mut_IDs = mut_IDs + ";" + str(j)
        mut_array.append([edgeID, mut_IDs])
        prev = prev + mut_num

    print_f(mut_array, out_f)
    return mut_array


def read_leaf_num(tree_f):
    ret = 0
    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    while(line != ""):
        line_a = line.split('\t')
        if line_a[5] == "1":
            ret = ret + 1
        line = file.readline().rstrip('\n')

    file.close()
    return ret 


def distribute_SNVcells(cell_n, tree_f, sigma, out_f):
    # to make all cell number add up to cell_n
    leaf_num = read_leaf_num(tree_f)
    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    SNVcell_array = []
    # SNV cell ID starts from 0 (0-based)
    prev = -1

    # to avoid negative number of cell, go through the whole file once and if negative, go through each line again
    leafnode = []
    while(line != ""):
        line_a = line.split('\t')
        # only look at leaves' percentages
        if line_a[5] == "1": 
            CNAcellP = float(line_a[4])
            leafnode.append([line_a[2], CNAcellP])
        line = file.readline().rstrip('\n')

    file.close()

    tag = True
    SNVcell_nums = init_m(len(leafnode), 2, 0)
    while tag:
        prev = 0
        for it in range(len(leafnode)):
            leafID, CNAcellP = leafnode[it]
            SNVcell_num = round(cell_n * (CNAcellP + random.gauss(0, sigma)))

            if SNVcell_num < 0:
                SNVcell_num = 0

            # 06132022 Add a protection so that the total will not be above cell_n
            if prev + SNVcell_num > cell_n or it == len(leafnode) - 1:
                SNVcell_num = cell_n - prev

            #print("In it " + str(it) + ", SNVcell_num = " + str(SNVcell_num) + ", total = " + str(prev + SNVcell_num))

            SNVcell_nums[it] = [leafID, SNVcell_num]
            if SNVcell_num < 0:
                print("In distributing SNV cells, SNVcell_num becomes negative. Start again. ")
                break
            prev = prev + SNVcell_num
            if it == len(leafnode) - 1:
                print("Success in distributing SNV cells. ")
                tag = False

    prev = -1
    for i in range(len(SNVcell_nums)):
        leafID, SNVcell_num = SNVcell_nums[i]
        SNVcell_IDs = "NA"
        if SNVcell_num != 0:
            SNVcell_IDs = str(prev + 1)
            for j in range(prev + 2, prev + SNVcell_num + 1):
                SNVcell_IDs = SNVcell_IDs + ";" + str(j)
        SNVcell_array.append([leafID, SNVcell_IDs])
        prev = prev + SNVcell_num


    print_f(SNVcell_array, out_f)
    return SNVcell_array


def add_missing(G, n, m, missingP):
    # in G (n*m) matrix, select missingP * n * m entries and flip them to 3. 
    missing_arr = random.sample(range(n * m), int(n * m * missingP))
    for i in missing_arr:
        row = int(i / m)
        column = i % m
        G[row][column] = 3
    return G


def count_total_value(M, n, m, value):
    total = 0
    for i in range(n):
        for j in range(m):
            if M[i][j] == value:
                total = total + 1
    return total


def add_FPFNs(D, n, m, alpha, beta, total_zeros, total_ones):
    FP_arr = random.sample(range(total_zeros), int(total_zeros * alpha))
    FN_arr = random.sample(range(total_ones), int(total_ones * beta))
    D_ = copy.deepcopy(D)
    index_neg = 0
    index_pos = 0
    for i in range(n):
        for j in range(m):
            if D[i][j] == 0:
                if index_neg in FP_arr:
                    # flip it
                    D_[i][j] = 1
                index_neg = index_neg + 1
            elif D[i][j] == 1:
                if index_pos in FN_arr:
                    # flip it
                    D_[i][j] = 0
                index_pos = index_pos + 1
    return D_
    


def reveal_edge_to_SNVcells(SNVcell_array, e_dict, n_dict, edge_revealP, edge_reveal_file, n):
    edge_reveal_f = open(edge_reveal_file, "w")
    # the edges selected
    edge_array = sorted(list(e_dict.keys()))
    # total edges to be revealed
    revealN = int(edge_revealP * len(edge_array))
    count = 0
    revealCells = []
    while count < revealN:
        edges = random.sample(range(len(edge_array)), 1)
        edge = str(edges[0])
        # find all cells below this edge separated by semicolon
        leaves = leaves_under_edge(edge, e_dict, n_dict)
        # find all the cell IDs on these leaves
        cellsIDs = []
        for j in range(len(SNVcell_array)):
            if SNVcell_array[j][0] in leaves:
                if SNVcell_array[j][1] != "NA":
                    # all those cell IDs on this row of SNVcell_array shall be recorded
                    for k in SNVcell_array[j][1].split(";"):
                        cellsIDs.append(k)
        if len(cellsIDs) == 0 or len(cellsIDs) > int(0.5 * n) or cellsIDs in revealCells:
            continue
        else:
            print("\t".join([edge, ";".join(cellsIDs)]), file=edge_reveal_f)
            count += 1
            revealCells.append(cellsIDs)
    edge_reveal_f.close()




# from tree_f we know which cells (on leaf nodes) have which mutations (on edges)
# n_ and SNVcell_array_ are for SCARLET. Not meaningful if forSCARLET is false.
def mutation_matrix(mut_array, SNVcell_array, tree_f, n, m, lossMutP, missingP, alpha, beta, edge_revealP, mutLoss_file, mutCNAoverlap_file, G_matrix_file, D_matrix_file, edge_reveal_file, D_miss_file, forSCARLET, n_, SNVcell_array_):

    # read mut_array to a dict whose first column is the edge IDs and second are the mutation IDs 
    mut_array_dict = {}
    for i in range(len(mut_array)):
        mut_array_dict[mut_array[i][0]] = mut_array[i][1]

    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    # for each leaf node, record all the edges that are above it
    # for each edge, record all the cells are below it


    # node_dict: a dict that has the child node ID as the key, .e and .p as the values, the edge above it and the parent node ID, respectively
    # edge_dict: a dict that has the edge ID as the key, .c and .p as the values, the child and parent node IDs on the two ends, respectivley
    edge_dict = {}
    node_dict = {}
    # c is a temporary data structure for retrieving all children nodes for all parent nodes for node_dict
    c_dict = {}
    while(line != ""):
        line_a = line.split('\t')
        # line_a[1] is parent, line_a[2] is child
        edge_dict[line_a[0]] = Edge(line_a[1], line_a[2])
        node_dict[line_a[2]] = Node(line_a[0], line_a[1], "NA")
        # record the children for each parent to be added to node_dict
        if line_a[1] not in c_dict.keys():
            c_dict[line_a[1]] = line_a[2]
        else:
            c_dict[line_a[1]] = c_dict[line_a[1]] + ";" + line_a[2]
        line = file.readline().rstrip('\n')

    file.close()

    for i in node_dict.keys():
        if i in c_dict.keys():
            node_dict[i].c = c_dict[i]
        else:
            # this is a leaf
            node_dict[i].c = "NA"
 

    # cell to mutation ID arrays
    SNVcell_mut = cell2mutation_arrays(SNVcell_array, node_dict, mut_array_dict) 

    SNVcell_mut_ = {}
    if forSCARLET:
        SNVcell_mut_ = cell2mutation_arrays(SNVcell_array_, node_dict, mut_array_dict) 

    # construct G matrix. n*m: n(cell), m(mut)
    G = makeG(n, m, SNVcell_mut)

    total_zeros = count_total_value(G, n, m, 0)
    total_ones = count_total_value(G, n, m, 1)
    #print("Before adding noises, total zeros for G: " + str(total_zeros) + "; total ones for G: " + str(total_ones))

    G_ = init(n_, m)
    if forSCARLET:
        G_ = makeG(n_, m, SNVcell_mut_)

    # loss of mutation, select lossMutP of the edges, for each edge, select without return a mutation that will be lost on cells below this edge.  
    G, mutEdgePair = addLossMut(G, lossMutP, edge_dict, node_dict, mut_array, SNVcell_array, mutLoss_file, mutCNAoverlap_file, m)

    
    if forSCARLET:
        # modify G_ according to the loss of mutation in nonSCARLET case
        G_ = addLossMutSCARLET(G_, edge_dict, node_dict, SNVcell_array_, n,  m, G, mutEdgePair)
            
    print_matrix(G, n, G_matrix_file + ".csv")

    if forSCARLET:
        print_matrix(G_, len(G_), G_matrix_file + ".SCARLET.csv")

        total_zeros = count_total_value(G_, n_, m, 0)
        total_ones = count_total_value(G_, n_, m, 1)
        print("Before adding noises, total zeros for G_: " + str(total_zeros) + "; total ones for G_: " + str(total_ones))

    # Step a. now make missing data to produce D 
    D_miss = add_missing(G, n, m, missingP) 
    print_matrix(D_miss, n, D_miss_file + ".csv")
    total_threes = count_total_value(D_miss, n, m, 3)
    #print("Total three for D_miss is " + str(total_threes))

    if forSCARLET:
        # copy the D_miss data for the first n rows, then call the add_missing for the rest of the rows
        D_miss_ = add_missing_SCARLET(D_miss, G_, n, m, missingP)
        print_matrix(D_miss_, n_, D_miss_file + ".SCARLET.csv") 
        total_threes = count_total_value(D_miss_, n_, m, 3)
        print("Total three for D_miss_ is " + str(total_threes))


    # Step b, c. flip 0's by \alpha and 1's by \beta 
    total_zeros = count_total_value(D_miss, n, m, 0)
    total_ones = count_total_value(D_miss, n, m, 1)
    #print("total zeros after missing: " + str(total_zeros) + "; total ones after missing: " + str(total_ones))
    D_miss_FP_FN = add_FPFNs(D_miss, n, m, alpha, beta, total_zeros, total_ones)
    #total_zeros = count_total_value(D_miss_FP_FN, n, m, 0)
    #total_ones = count_total_value(D_miss_FP_FN, n, m, 1)
    #print("total zeros after FPFN: " + str(total_zeros) + "; total ones after FPFN: " + str(total_ones))

    # store the D matrix
    print_matrix(D_miss_FP_FN, n, D_matrix_file + ".csv")

    if forSCARLET:
        # copy the D_miss_FP_FN matrix for the first n rows, then call the add_FPFNs for the rest of the rows
        D_miss_FP_FN_ = add_FPFNs_SCARLET(D_miss_FP_FN, D_miss_, n, m, alpha, beta)
        print_matrix(D_miss_FP_FN_, n_, D_matrix_file + ".SCARLET.csv")
        total_zeros = count_total_value(D_miss_FP_FN_, n_, m, 0)
        total_ones = count_total_value(D_miss_FP_FN_, n_, m, 1)
        print("After adding noises, total zeros for G_: " + str(total_zeros) + "; total ones for G_: " + str(total_ones))

    if forSCARLET:
        print("Success in generating D and G matrix for SCARLET. ")
    
    # Lastly, reveal the edges to some SNV cells. 
    reveal_edge_to_SNVcells(SNVcell_array, edge_dict, node_dict, edge_revealP, edge_reveal_file, n)


   
# beginning of the program
if len(sys.argv) <= 1:
    print("""
    This generates the mutation matrix with the ground truth data. 
    Usage: python sim_par.py -a [alpha] -b [beta] -m [missing-rate] -l [loss-rate] -c [num-cells] -n [num-mutations] -p [perc-CNAedge-SNVcells] -s [sigma] -f [input-tree-file] -P [prefix-output-files]
        -a (--alpha)        False positive rate. [0.01]
        -b (--beta)         False negative rate. [0.2]
        -m (--missing-rate) Missing rate in G. [0.2]
        -l (--loss-rate)    Mutation loss rate due to copy number loss (rate is on edge number). [0.1]
        -c (--num-cells)    Total number of SNV cells. [100]
        -n (--num-mutations)    Total number of mutations. [100]
        -p (--perc-CNAedge-SNVcells)    Percentage of CNA edges revealed to SNV cells. [0.1]
        -s (--sigma)        Standard deviation between the SNV and CNA trees on the number of cells each leaf carries. [0.05]
        -f (--tree-file)    The input tree structure file. ["NA"]
        -P (--prefix)       Prefix of output files. 
        -S (--if-scarlet)   Turn this option on if SCARLET input will also be generated. Extra cells (one for each internal node) will be attached to the .G, .miss, .D matrices if this is turned on. [True]
    """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='This script iterates all parameters to generate simulated data.')
parser.add_argument('-a', '--alpha', default=0.01)
parser.add_argument('-b', '--beta', default=0.2)
parser.add_argument('-m', '--missing-rate', default=0.2)
parser.add_argument('-l', '--loss-rate', default=0.1)
parser.add_argument('-c', '--num-cells', default=100)
parser.add_argument('-n', '--num-mutations', default=100)
parser.add_argument('-p', '--perc-CNAedge-SNVcells', default=0.1)
parser.add_argument('-s', '--sigma', default=0.05)
parser.add_argument('-f', '--tree-file', default="NA")
parser.add_argument('-P', '--prefix', default="NA")
parser.add_argument('-S', '--if-scarlet', action='store_true')


args = parser.parse_args()
alpha = float(args.alpha)
beta = float(args.beta)
missing = float(args.missing_rate)
loss = float(args.loss_rate)
cell_n = int(args.num_cells)
mut_n = int(args.num_mutations)
edge_reveal = float(args.perc_CNAedge_SNVcells)
sigma = float(args.sigma)
tree_f = args.tree_file
prefix = args.prefix
forSCARLET = args.if_scarlet

# First, distribute the mutations on the edges. Read tree_f output to prefix + ".mut.csv" 
mut_array = distribute_mutations(mut_n, tree_f, prefix + ".mut.csv")

# Second, distribute the cells on the leaf nodes. Read tree_f output to prefix + ".SNVcell.csv" 
SNVcell_array = distribute_SNVcells(cell_n, tree_f, sigma, prefix + ".SNVcell.csv")

# for SCARLET input, we add one SNV cell to all internal nodes, and those leaves that do not have any nodes. 
SNVcell_array_SCARLET = []
total_leaves_SCARLET = 0
if forSCARLET:
    SNVcell_array_SCARLET, total_leaves_SCARLET = add_SNVcells_internal(tree_f, SNVcell_array, cell_n, prefix + ".SNVcell.SCARLET.csv")

# Last, make mutation matrices G and D. 
# .miss.csv, .G.csv and .D.csv for non-SCARLET input. .miss.SCARLET.csv, .G.SCARLET.csv and .D.SCARLET.csv for SCARLET input. No .reveal file for scarlet. The same .mutCNAoverlap.csv and .mutLoss.csv file for both non-SCARLET and SCARLET input. 
mutation_matrix(mut_array, SNVcell_array, tree_f, cell_n, mut_n, loss, missing, alpha, beta, edge_reveal, prefix + ".mutLoss.csv", prefix + ".mutCNAoverlap.csv", prefix + ".G", prefix + ".D", prefix + ".reveal.csv", prefix + ".miss", forSCARLET, total_leaves_SCARLET, SNVcell_array_SCARLET)



