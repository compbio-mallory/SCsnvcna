Note: This repository contains scripts that were developed for "SCsnvcna: Integrating SNVs and CNAs on a phylogenetic tree from single-cell DNA sequencing data". 

Authors: Liting Zhang, Hank Bass, Jerome Irianto, and Xian Mallory. 

# SCsnvcna
SCsnvcna is an algorithm that integrate SNVs and CNAs on a phylogenetic tree from single-cell DNA sequencing data. SCsnvcna uses a loss-supported model that allow mutation lossess based on observed copy number loss on CNAs tree. 


SCsnvcna is implemented in Python. 

## Contents
1. [Setup](#setup) 
	- [Dependencies](#dependencies)

2. [Running SCsnvcna](#runningscsnvcna)
	- [Input](#input)
	- [Output](#output)
	- [Usage](#usage)
3. [Simulation](#getsimulation)

<a name="setup"></a>
## Setup


<a name="dependencies"></a>
### Dependencies

- Python 3.9 ([anaconda distribution](https://www.anaconda.com/distribution/) recommended)  


### SCsnvcna Setup
To use SCsnvcna, first clone the SCsnvcna repository locally. 

```git clone git@github.com:compbio-mallory/SCsnvcna.git```

SCsnvcna uses Python and doesn't require any compilation. 


<a name="runningscsnvcna"></a>
## Running SCsnvcna
<a name="input"></a>
### Input Files
SCsnvcna takes as input three files. The first describe the copy number tree structure from CNAs cell. The second descibes the observed genotype matrix from SNV cells. The third describes the set of supported losses for each edge in the copy number tree. 

1. **Tree structure file**. This file describe the copy-number tree structure constructured from CNAs cells. This file takes the format of a tab-separated edge list, where the first two columns are child and parent nodes, followed by  edge ID connecting into child node (same as child node), edge length, percentage of CNAs cells below this node, whether it is a leaf, and an optional node name. 
```
 0	-1	0	0.001	1	0
 1	0	1	0.002	0.6	0
 2	0	2	0.04	0.4	1
 3	1	3	0.05	0.3	1
 4	1	4	0.4	0.3	1
```
2. **Genotyp matrix file**. This file describe the observed genotype matrix from SNVs cells. This file takes the format of a tab-seperated matrix. Where each line is the genotype for each cell, whereas each columns are genotype for each SNVs. 0 is absent, 1 is present, and 3 is missing.
 ```
 0	1	3	0	...
 1	1	1	0	...
 ...
 ```
3. **Mutation loss support file**. This file describe the potential mutation loss of each edge. The first column is edge, the second column is the SNV that might be lost on this edge. 
 ```
 1 19
 12 55
 ...
 ```	
 
Examples of these files can be found in the `data/simulation/` directory.

<a name="output"></a>
### Output Files

Scsnvcna produces several output files, as detailed below. In addition to these, SCsnvcna produces several auxiliary files that are used for graphing purposes and log file are not detailed below. 

1. **Binary mutation matrix file** (`[prefix].G`). This file describes the binary presence (`1`) or absence (`0`) of each mutation in each cell. This is a tab-separated file where rows correspond to cells, of the following format. 
	```
	[mut1] [mut2] ... [mut_m]
	...
	```
2. **SNV placement file** (`prefix.snv`). This file describe the placement of SNVs. The first column is the SNV ID, and the second column is where the SNV is placed on.
	```
	SNV1	Edge1
	...
	```
3. **SNV cell placement file** (`prefix.cell`). This file describe the placement of SNVs cell. The first column is the SNV cell id, and the second column is where the SNV cell is placed on.
 	```
	cell1	leave1
	...
	```
	
4. **Log probability of the selected tree** (`prefix.Prob`). This file describe probability of the selected tree.
 	```
	Probability
	```
### Usage

SCsnvcna can be run from the command line as follows.

```
python3 code/main.py -tree data/sample_data/tree_8_p3.tsv -D data/sample_data/input.D.tsv -overlap data/sample_data/input.mutCNAoverlap.tsv -out data/sample_data/prefix -alpha 0.01 -beta 0.2 -sigma 0.05 -searches 100000 -reveal data/sample_data/input.SNVcell.tsv -restart 10
```
 - -tree Copy number tree structure file
 - -D observed genotype matrix for SNV cells
 - -overlap files contain potential mutation loss. Empty file if no mutation loss. 
 - -out output prefix 
 - -alpha initial false positve rate
 - -beta initial false negative rate
 - -sigma standard deviation of CPs between CNA cells and SNV cells
 - -searches number of iteration. Default is 100000. 
 - -reveal file contains where should SNV cells be placed. Empty file if no SNV cells constrian. 
 - -restart number of restart. Default is 10.

To visualize the placement of SNVs and SNV cells
```
Rscript code/plot_tree.R prefix.treeText prefix
```
 - prefix.treeText are the newick output of SCsnvcna.
 - prefix output plot prefix

To tune the plot, change the parameter at then end of the `plot_tree.R` accordingly.  

<a name="getsimulation"></a>
## Get Simulation Data for SCsnvcna
1. Generate a Tree
```
python3 code/simulation_code/gen_tree.py -F 8 -B 0.3 -o tree_8_p3.tsv
```
 - -F width of the tree. Number of leaves on this tree. Default is 16.
 - -B Beta splitting variable. Default is 0.3.
 - -o Output file. 
2. Add mutations and cells to the tree.
```
python3 code/simulation_code/sim_par.py -c 50 -n 20 -f tree_8_p3.csv -P input
```
 - -c Number of cells. Default is 100.
 - -n Number of mutation. Default is 100.
 - -f Tree file from step 1.
 - -P Output file prefix. 
