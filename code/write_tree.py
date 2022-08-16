from cProfile import label
import os, sys
from numpy import roots
from numpy.core.fromnumeric import sort
from numpy.lib.function_base import append
from numpy.lib.shape_base import split
from data_structure import MySNV1, MyNode1
from random import sample

import numpy as np

# function to write out tree information
def writeTree(f, tree):

    G = tree.G
    for row in G:
        f.write(' '.join([str(a) for a in row]) + "\n")
    f.write("SNV placement\n")
    for j in tree.snvs:
      f.write(str(j.id) + "," + str(j.edge) + "," + str(j.loss) + "\t")
    f.write("\n")
    


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
  snvs = ''
  for snv in nodes[root].new_snvs:
    snvs = snvs + str(snv) + "|"
  rootNHX = ''
  if len(snvs) > 1:
    rootNHX = "):" + str(round(nodes[root].edge_length, 2)) +"[&&NHX:S=" + snvs + ":N=" + str(rootName) + "];"
  else:
    rootNHX = "):" + str(round(nodes[root].edge_length, 2)) + "[&&NHX:" + ":N=" + str(rootName) + "];"
  numChild = len(nodes[root].children)
  for c in nodes[root].children:
    #print(root)
    #print(c)
    child = nodes[c]
    #print(isinstance(nodes[c].name, int))
    rootNHX = genTreeTextChild(child, nodes) + rootNHX
    numChild -= 1
    if(numChild > 0):
      rootNHX = "," + rootNHX
  rootNHX = "(" + rootNHX +"\n"
  return(rootNHX)

def genTreeTextChild(node, nodes):
  snvs = '+'
  for snv in node.new_snvs:
    snvs = snvs + str(snv) + "|"
  snvsLoss = "-"
  if node.loss_snvs:
    #print("mutation loss")
    for snv in node.loss_snvs:
      snvsLoss = snvsLoss + str(snv) + "|"
  cells = ''
  for cell in node.cells:
    cells = cells + str(cell) + "|"
  currText = ''
  if(node.if_leaf):
    currText = cells + ":"
  e_len = round(node.edge_length + 0.1, ndigits= 2)

  if node.edge_length == 0:
    e_len = 0
  currText =  currText  + str(e_len) + "[&&NHX:"
  if(len(snvs) > 1 or len(snvsLoss) > 1):
    currText = currText + ":S="
    if(len(snvs) > 1):
      currText = currText + snvs
    if(len(snvsLoss) > 1):
      currText = currText  + snvsLoss

  cna = "|"
  counter = 0  
  gain = 0
  loss = 0 
  
  node_label = str(node.name)
  # for raw paup tree
  # if isinstance(node.name, int):
  #   node_label = str(int(node.name) + 2)    
    
  gain_ = []
  loss_ = []
  for i in range(len(node.cn_summary)):
    for key in node.cn_summary[i]:
      s = int(key.split(".")[0])
      end = int(key.split(".")[1])
      curr_cn = node.cn_summary[i][key] 
      prev = node.parentID 
      prev_cn = 2
      if(key in nodes[prev].cn_summary[i]):
        prev_cn = nodes[prev].cn_summary[i][key]
      if(curr_cn != prev_cn):
        if(curr_cn - prev_cn > 0):
       
          cna = cna + "+" + str(counter) + "|"
          gain += (end - s)
          gain_.append([counter, i+1, s, end])
        if(curr_cn - prev_cn < 0):
     
          cna = cna + "-" + str(counter) + "|"
          loss += (end - s)
          loss_.append([counter, i+1, s, end])
      counter += 1

  currText = currText #+ ":M=" + cna #+ ":G=" + str(gain) + ":L=" + str(loss)
  #currText = currText +  ":N=" + str(node.name) + "]"
  
  # node_label = str(node.id)
  currText = currText +  ":N=" + (node_label) + "]"
  if(len(node.children) > 0):
    currText = "):" + currText
    numChild = len(node.children)
    #print("children ", node.children)
    for c in node.children:
      #print("get node ", c)
      child = nodes[c]
      text = genTreeTextChild(child, nodes)
      currText = text + currText
      numChild -= 1
      if(numChild > 0):
        currText = "," + currText
    currText = "(" + currText
  return currText




def getError(D, G):
    fn = 0
    fp = 0
    num1 = 0
    num0 = 0
    for i in range(len(D)):
      for j in range(len(D[0])):
        if(D[i][j] == 1):
          if(G[i][j] == 0):
            fp += 1
            num0 += 1
          else:
            num1 +=1
        if(D[i][j] == 0):
          if(G[i][j] == 1):
            fn +=1
            num1 += 1
          else:
            num0 += 1
    FP = 0
    FN = 0
    if num0 != 0:
      FP = fp / num0
    if num1 != 0:
      FN = fn / num1
    return FP, FN


