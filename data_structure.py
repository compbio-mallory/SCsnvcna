# definition of one SNV
class MySNV1():
    def __init__(self, id, name, ale, chr, pos, ori_nuc, new_nuc):
        #MySNV.__init__(self, ale, chr, pos, ori_nuc, new_nuc)
        self.id = id
        self.name = name
        self.ale = ale
        self.chr = chr
        self.pos = pos
        self.ori_nuc = ori_nuc
        self.new_nuc = new_nuc
        self.edge = -1 #where it is placed on 
        self.overlapped = []

class CN1:
    def __init__(self, CN_Ale, CN_Del, CN_chromosome, CN_p1, CN_p2, CN_amp_num, corres):
        self.CN_Ale = CN_Ale
        self.CN_Del = CN_Del #if deletion 1
        self.CN_chromosome = CN_chromosome
        self.CN_p1 = CN_p1
        self.CN_p2 = CN_p2
        self.CN_amp_num = CN_amp_num
        self.corres = corres
    def get_CN_Ale(self):
        return self.CN_Ale
    def get_CN_Del(self):
        return self.CN_Del
    def get_CN_chromosome(self):
        return self.CN_chromosome
    def get_CN_position(self):
        #print self.CN_p1, self.CN_p2
        return self.CN_p1, self.CN_p2
    def get_CN_amp_num(self):
        return self.CN_amp_num


# tree is an array of MyNode
class MyNode1:
    def __init__(self, id = -1, name = None, parent=None):
        self.id = id
        self.name = name
        self.parent = parent
        self.children = []# ids for children
    #    self.is_dead=False 
        self.edge_length = 1
    #    alelle length for each chromosome, root has the same as reference
    #    self.cn=[]#array of CN
    #    self.chrlen=[]#chromosome length
    #    self.ref=[]#ref[i][j], for each allele, each chromosome 
        self.snvs = []#snv ids 
    #    self.corres = []#corres[i][j][k], for each allele, each chromosome, each cnv region ref(r1, r2), gen(g1,g2)
        self.cn_summary = [{}] * 24 #for each chromosome, cn_summary[i] = {"s:e": cn, ...}, 
    #    self.cn_detail=[]#cn_detail[i][j], for each allele, each chromosome, 
        self.parentID = -1
    #    self.true_CNs = []
    #    self.depth_ = -1#?
        self.perc = -1#percentage of ícells below this node
        self.if_leaf = False
        self.cells = []
        self.new_snv = [] # new snvs happen on this node, snv ids
        self.loss_snv = []
    #    self.aberrations = []#?
    def getTuple(self):
        return self.tuple
    # def setDead(self):
    #     self.is_dead=True
    def getID(self):
        return self.id
    # def getDepth(self):
    #     return self.depth_
    def getPerc(self):
        return self.perc


# tree structure
class TREE():
    def __init__(self, id = -1, p = 0, nodes = [], G = [], D_bar = [], D_hat = [], theta = [], sigma = 0, snvs = [], cells = []):
        self.id = id
        self.p = p
        self.G = G
        self.nodes = nodes 
        self.D_bar = D_bar
        self.D_hat = D_hat
        self.theta = theta
        self.sigma = sigma
        self.snvs = snvs
        self.cells = cells

def get_Tree_P(tree):
    #print(tree.p)
    return tree.p[0]