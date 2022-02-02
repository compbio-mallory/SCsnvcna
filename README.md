# SCsnvcna
## Files
* `SCsnvcna.py` is the main program for SCsnvcna.</br>
* `data_structure.py` is the file containing basic data structure used in main program.
* `read_files.py` is the file containing functions to read input files.
* `CN_call.R` is the program to generate copy number profile.</br>
* `merge_seg_qdnaseq.py` is the program to merge all the copy number files. 
* `Nexus_gen.py` is the program to generate nexus file.
* `CN.py` is a necessary file to read the simulated data.
* `gen_tree_overlappingCNA.py` is a necessary file to read the simulated data. 
* `CN.py` is a necesary file to read the simulated data.
* `Gen_Ref_Fa.py` is a necessary file to read the simulated data. 
* `plot_tree.R` is a program to plot the resulting tree.

## Unify the CNA Tree and SNV
### Working with simulation data
1. Generate simulated single cell phylogenetic tree with copy number aberration and single nucletide variants on it.
2. `python3.9 SCsnvcna.py  --simulated 1 --treeFile from_first_step.tree.npy --snvs_table gt.snv.csv --alpha 0.0174 --beta 0.1256 --sigma 0.01 --imputedFP 0.01 --imputedFN 0.3 --imputedMiss 0.1 --out test1_01_3_1 --errorSD 0.1 --searches 100000`
  * `--simulated` set to 1 if working on simulation data. Default is 0.
  * `--treeFile` file that containing the CNA tree information. This file comes from the simulator program. 
  * `--snvs_table` file that containing the snv information per single cell. This file coms from the simulator program.
  * `--alpha` false positive rate. The default is 0.0174.
  * `--beta` false negative rate. The default is 0.1256.
  * `--sigma` the variance between cells generate the CNA tree and cells generate SNV profiels. The default is 0.001.
  * `--imputedFP` false positive rate that is used to impute the SNV table.
  * `--imputedFN` false negative rate that is used to impute the SNV table.
  * `--imputedMiss` missing rate that is used to impute the SNV table.
  * `--out` output file name.
  * `--errorSD` standard deviation that is used to impute an incorrect CNA tree. If `--errorSD` is set to be greater than 0, then the CNA tree will be altered based on `--errorSD`. Default is 0.
  * `--searches` the number of iteration. Default is 1000000.
  * The program will return three files `$out`, `$out_evaluation`, `$out_treeText`. `$out` contains log inforamtion for selected top 10 trees. `$out_evaluation` contains evaluation matrix for selected top 10 trees. `$out_treeText` contains text connection for each selected tree for plotting purpose.

### Working with real data

#### Generate Copy Number Profile and CNA Tree
1. Get segments data for each bam file using QDNAseq package using a bin width of 100kbp.`CN_call.R`. This step will return a tsv file containing the copy number profile for each cell and segment. 
2. Merge segments file from previous step. `merge_seg_qdnaseq.py $input $output`. This step find the common breakpoints and merges unique CNA segments. `$input` is the file from previous step. `$output` is the output file name. This program will return two files `$output` and `$output.unique`, which only contains segments that have copy number variation among cells.
3. Create Nexus file for paup. `Nexus_gen.py $input $output.nex`, where `$input` is the input file from previous step.
4. Use command `execute $nexFile` inside the paup interative program to generate a logfile.

#### Generate SNV Profile and SNV Table.
* SNV call can be done using samtools, GATK, freebayes or similar program.
* SNV profile should be a list of SNVS. See the sample SNV profile below:
  * LRP1B chr2  140851729
  * SNV chromosome  position
  * ... ... ...
* SNV table should looks like this
  * cell1 1 0 -1  ... 1
  * cell2 1 1 0 ... 0
  * cellX SNV1  SNV2  ...
  * ... ... ... ... ...
  * where 1 is present, 0 is absent, -1 is missing. 

#### Run SCsnvcna
`python3.9 SCsnvcna.py --cnFile merged_pa_100kbp.tsv.unique --treeFile merged_pa_100kbp.log --snvs mutations38filtered.txt --snvs_table snv_table_filtered38 --out pa_100kbp_0174_1256_001_1 --alpha 0.0174 --beta 0.1256 --sigma 0.001 --treeNum 1`
- $cnFile is the *.unique file from generate copy number profile step.
- $treeFile is the *.log file from generate copy number profile step.
- $snvs is the snvs profile file.
- $snvs_table is the snvs table for each cell.
- $outputFile is the ouptufile name.
- $initialAlpha is the initial false positive rate. Default is 0.0174. 
- $initialBeta is the initial false negative rate. Default is 0.1267.
- $initialSigma is the initial standard deviation for CNA tree and SNVs. Default is 0.001.
- $treeNumber is the tree number that should be selected from $treeFile to used. Default is 1. 
- The program will return three files `$out`, `$out_evaluation`, `$out_treeText`. `$out` contains log inforamtion for selected top 10 trees. `$out_evaluation` contains evaluation matrix for selected top 10 trees. `$out_treeText` contains text connection for each selected tree for plotting purpose.
