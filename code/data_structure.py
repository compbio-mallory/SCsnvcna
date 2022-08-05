
# SNV 
# definition of one SNV
class MySNV1():
  def __init__(self, id, name = -1, ale = -1, chr = -1, pos = -1):
    #MySNV.__init__(self, ale, chr, pos, ori_nuc, new_nuc)
    self.id = id
    self.name = name
    self.ale = ale
    self.chr = chr
    self.edge = -1 #where it is placed on 
    self.overlapped = []
    self.loss_node = [] #where this SNV will be lost
    self.loss = -1
    self.pos = pos
# tree is an array of MyNode
class MyNode1:
  def __init__(self, id = -1, name = None, parent=None, parentID = -1, perc = None, if_leaf = False):
    self.id = id
    self.name = name
    self.parent = parent
    self.children = []# ids for children
    self.edge_length = 1
    self.snvs = []#snv ids 
    self.cn_summary = [{}] * 24 #for each chromosome, cn_summary[i] = {"s:e": cn, ...}, 
    self.parentID = parentID
    self.perc = perc#percentage of ícells below this node
    self.if_leaf = if_leaf
    self.cells = []
    self.new_snvs = [] # new snvs happen on this node, snv ids
    self.loss_snvs = None


