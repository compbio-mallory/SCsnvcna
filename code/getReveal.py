import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import umap
import hdbscan 
from collections import Counter

SNVcell = sys.argv[1]# snv cell list with id
CNAcell = sys.argv[2] # input CNA tree with id and names for CNA cells
SNVCN = sys.argv[3]
CNACN = sys.argv[4]

# get bins that overlap CNAs cells and SNVs cells
def getBins(SNVdf, CNAdf):
    # get segments that overlaps for 90%
  overlap_dict = {}
  for index_S, row_S in SNVdf.iterrows():
    CHROM = row_S["CHROM"]
    START_S = row_S["START"]
    END_S = row_S["END"]
    # Filter dataframe CNA by CHROMosome
    C_CHROM = CNAdf[CNAdf["CHROM"] == CHROM]
    # Calculate the length of the segment in A
    segment_length_S = END_S - START_S
    # Find the rows in dataframe B that overlap for more than 80%
    for index_C, row_C in C_CHROM.iterrows():
      START_C = row_C["START"]
      END_C = row_C["END"]
      if END_C < START_S:
        continue # current segment in CNA is before SNV, next
      if START_C > END_S:
        break # no need to move forward
      # Calculate the length of the overlap between SNV and CNA
      overlap_length = min(END_S, END_C) - max(START_S, START_C)
      # Calculate the length of the segment in C
      segment_length_C = END_C - START_C
      # Check if the overlap is more than 80% of the segment in A and B
      if overlap_length >= 0.9 * segment_length_S and overlap_length >= 0.9 * segment_length_C:
        overlap_dict[index_S] = index_C
        #print(SNVdf.iloc[index_S]["START"], SNVdf.iloc[index_S]["END"] )
        SNVdf["START"].iloc[index_S] = CNAdf["START"][index_C]
        SNVdf["END"].iloc[index_S] = CNAdf["END"][index_C]
        #print(SNVdf.iloc[index_S]["START"], SNVdf.iloc[index_S]["END"] )
  return SNVdf

def group_rows_with_gaps(df):
  result_df = pd.DataFrame(columns=['CHROM', 'START', 'END'] + list(df.columns[3:]))
  current_group = None
  for index, row in df.iterrows():
    #print(current_group)
    if current_group is None:
      current_group = row
    else:
      prev_END = current_group['END']
      if prev_END + 1== row['START']:
        current_group['END'] = row['END']
        for column in df.columns[3:]:
          if row[column] != current_group[column]:              
            result_df = result_df.append(current_group, ignore_index=True)
            #result_df = pd.concat([result_df, current_group], axis=0)
            current_group = row
            break
      else:
        result_df = result_df.append(current_group, ignore_index=True)
        #result_df = pd.concat([result_df, current_group], axis=0)
        current_group = row
  #result_df = pd.concat([result_df, current_group], axis=0)
  result_df = result_df.append(current_group, ignore_index=True)
  return result_df

def find_most_common_values(lst):
	# Get count of each element in the list
	count_dict = Counter(lst)

	# Get the most common value and its count
	most_common = count_dict.most_common(1)[0]
	most_common_value = most_common[0]
	most_common_count = most_common[1]

	# Calculate the percentage of times the most common value occurs in the list
	perc = (most_common_count / len(lst)) * 100
	if perc >= 90:
		return most_common_value, most_common_count
	else:
		return -1, -1
        
# find columns specific CNAs
def getCloneCNAs(CNAdf, clusters, SNVdf):
  print(CNAdf.head())
  print(SNVdf.head())
  merged_df = pd.merge(CNAdf, SNVdf, on=['CHROM', 'START', 'END'])
  CNAcellsCols = CNAdf.columns[:-3].tolist()
  SNVonly = merged_df.drop(CNAcellsCols, axis = 1)
  CNAonly = merged_df.loc[:, CNAcellsCols+ ['CHROM', 'START', 'END']]
  uniqueCNAs = {}
  for c in clusters:
    df = pd.DataFrame(columns=['CHROM', 'START', 'END', 'CN'])
    uniqueCNAs[c] = df
  for index, row in merged_df.iterrows():
    for c in clusters:
      cells = clusters[c]
      cns = row[cells].tolist()
      val, count = find_most_common_values(cns)
      if val == -1 or val == 2:
        continue
      # how many other cells have this CNA
      num_freq = CNAonly.iloc[index].value_counts()[val] - count
      if num_freq < 0.05 *(CNAonly.shape[1] - 3 - len(cells)):
        # unique CNA to this group
        #print(c, index, num_freq, val,CNAonly.shape[1] - 3 - len(cells))
        newRow = {"CHROM":row["CHROM"], "START":row['START'],"END": row['END'], "CN":val}
        uniqueCNAs[c] = uniqueCNAs[c].append(newRow, ignore_index=True)
        
  SNV2CNA = {} #how many CNAs does a SNV cell contain for each CNA cluster
  for c in uniqueCNAs:
    df = pd.merge(uniqueCNAs[c], SNVonly, on =['CHROM', 'START', 'END'])
    #print(df.head)
    for index, row in df.iterrows():
      cn_value = df.at[index, "CN"]
      same_cn_cols = df.columns[df.iloc[index] == cn_value].tolist()
      for col in same_cn_cols:
        if col not in SNV2CNA:
          SNV2CNA[col] = [0] * len(clusters)
          SNV2CNA[col][c] = 1
        else:
          SNV2CNA[col][c] += 1

  return uniqueCNAs, SNV2CNA

SNVdf = pd.read_csv(SNVCN, sep="\t")
CNAdf = pd.read_csv(CNACN, sep="\t")

SNVdfCP = getBins(SNVdf, CNAdf)
# remove common CNAs from CNAs cells to reduce run time
CNAdfCP = CNAdf.drop(["CHROM","START","END"], axis = 1)

def check_row(row):
    return row.value_counts(normalize=True).iloc[0] > 0.6
# Apply the function to each row of the dataframe
CNAdfCP["filter"] = CNAdf.apply(check_row, axis=1)
CNAdfCP["CHROM"] = CNAdf["CHROM"]
CNAdfCP["START"] = CNAdf["START"]
CNAdfCP["END"] = CNAdf["END"]
CNAdfCP = CNAdfCP[CNAdfCP['filter'] != True]
CNAdfCP = CNAdfCP.drop(["filter"], axis = 1)
print(CNAdfCP.head())
# cluster CNA cells
CN = CNAdfCP.drop(["CHROM","START","END"], axis = 1)
CN = CN.values.T
clusterable_embedding = umap.UMAP(
		n_neighbors=3,
		min_dist=0.0,
		n_components=3,
		random_state=42,
).fit_transform(CN)
labels = hdbscan.HDBSCAN(
	min_samples=3,
	min_cluster_size=3,
	metric = "manhattan"
).fit_predict(clusterable_embedding)

# get cluster for each cell
clusters = {}
CNAcells = CNAdf.columns.to_list()[3:]
CNAcells = np.array(CNAcells)
for i in range(max(labels) + 1):
    ind = np.where(labels == i)
    clusters[i] = CNAcells[ind[0]]

uniqueCNAs, SNV2CNA = getCloneCNAs(CNAdfCP, clusters, SNVdfCP)
revealDict = {}
for c in SNV2CNA:
  if c == "CN":
    continue
  for i in range(len(SNV2CNA[c])):
    print(c, SNV2CNA[c][i], len(uniqueCNAs[i]) )
    if len(uniqueCNAs[i]) < 10:
      continue
    if SNV2CNA[c][i]/ len(uniqueCNAs[i]) >= 0.8:
      if c not in revealDict:
        revealDict[c] = clusters[i].tolist()
      else:
        revealDict[c].extend(clusters[i].tolist())

f = open(CNAcell, "r")
CNAcell2id = {}
lines = f.readlines()
CNAcell2id = {}
for line in lines:
	cell = line.rstrip().split()[6]
	id = line.rstrip().split()[0]
	CNAcell2id[cell] = id
f.close()

f = open(SNVcell,"r")
lines = f.readlines()
SNVcell2id = {}
for line in lines:
	cell = line.rstrip().split()[1]
	id = line.rstrip().split()[0]
	SNVcell2id[cell] = id
f.close()

#f = open("reveal1.tsv", "w")
for key in revealDict:
	SNVid = SNVcell2id[key]
	for v in revealDict[key]:
		if v not in CNAcell2id:
			continue
		CNAid = CNAcell2id[v]
		print(str(CNAid) + "\t" + str(SNVid) )
		#f.write(str(CNAid) + "\t" + str(SNVid) + "\n")
#f.close()