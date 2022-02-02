import sys, os, math, re

# take a split line from input copynumber file
# transfer copy number back to integer
# @segs a line from input copy number
# return a list [chrom, start, end, cn1, cn2 ...]
def getCopyNumber(line):
    line = line.rstrip().split()
    if(line[1] == "X"):
        line[1] = "23"
    if(line[1] == "Y"):
        line[1] = "24"
    res = [int(line[1]), int(line[2]), int(line[3])]
    for i in range(4, len(line)):
        cn = round(2**float(line[i])*2)
        res.append(cn)
    return res

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    f = open(infile, "r")
    OUT = open(outfile, "w")
    header = f.readline().split()
    OUT.write("CHROM\tSTART\tEND\t")
    for i in range(4, len(header)):
        sample = re.sub('[^A-Za-z0-9]+', '', header[i])
        OUT.write(sample + "\t")
    OUT.write("\n")
    firstLine = f.readline() 
    prev = getCopyNumber(firstLine)
    start = prev[1]
    end = prev[2]
    curr = []
    segs = []
    for line in f:
        curr = getCopyNumber(line)
        if(curr[0] != prev[0]):# if different chromosome 
            OUT.write(str(prev[0]) + "\t")
            OUT.write(str(start) + "\t")
            OUT.write(str(end) + "\t")
            temp = [prev[0], start, end]
            for i in prev[3:]:
                OUT.write(str(i) + "\t")
            OUT.write("\n")
            start = curr[1]
            end = curr[2]
            prev = curr
        elif(curr[1] - prev[2] > 500001):#if not continue segs
            OUT.write(str(prev[0]) + "\t")
            OUT.write(str(start) + "\t")
            OUT.write(str(end) + "\t")
            for i in prev[3:]:
                OUT.write(str(i) + "\t")
            OUT.write("\n")
            start = curr[1]
            end = curr[2]
            prev = curr
        else:
            for i in range(3, len(prev)): #if read number is different
                if(curr[i] != prev[i]):
                    OUT.write(str(prev[0]) + "\t")
                    OUT.write(str(start) + "\t")
                    OUT.write(str(end) + "\t")
                    for i in prev[3:]:
                        OUT.write(str(i) + "\t")
                    OUT.write("\n")
                    start = curr[1]
                    end = curr[2]
                break
            end = curr[2] # continue and same cn for all cells, merge seg
            prev = curr

    OUT.write(str(prev[0]) + "\t")
    OUT.write(str(prev[1]) + "\t")
    OUT.write(str(prev[2]) + "\t")
    for i in prev[3:]:
        OUT.write(str(i) + "\t")
    OUT.write("\n")
    OUT.close()
    f.close()
    outfile1 = outfile + ".unique"
    f = open(outfile, "r")
    OUT = open(outfile1, "w")
    header = f.readline()
    OUT.write(header)
    for line in f:
        temp = line.rstrip().split()
        cur = temp[3]
        for cn in temp[3:]:
            if cn != cur:
                OUT.write(str(temp[0]) + "\t" + str(temp[1]) + "\t" + str(temp[2]) + "\t")
                for cn1 in temp[3:]:
                    OUT.write(str(cn1) + "\t")
                OUT.write("\n")
                break
if __name__ == "__main__":
  main()