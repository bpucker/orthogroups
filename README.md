# orthogroups
Scripts related to the filtering and analysis of orthogroups in the Caryophyllales. [OrthoFinder2](https://github.com/davidemms/OrthoFinder) was used for the initial identification of orthogroups.


### 1) Remove orthogropus with low taxa represenation (likely contamination/artifacts) ###

This script takes the orthogroup.tsv file produced by OrthoFinder2 and transfers only orthogroup with a certain level of taxa representation into a new file. The --optimize options evaluates different cutoff values and supports the identification of a suitable one.


```
python filter_OG_by_taxa.py
--in <INPUT_OG_FILE>
--out <OUTPUT_OG_FILE>
--taxa <MINIMAL_TAXA_NUMBER> | --optimize
```          

`--in` specifies the OrthoFinder2 result while which contains the orthogroup sequence IDs (TSV).

`--out` specifies the output file of this filtering process.

`--taxa` specifies the minimal number of taxa to be reprented in an orthogroup to keep the group. Orthogroups with less than this number of taxa are dropped.

`--optimize` activates an iterative process which performs the filtering for different taxa cutoffs and outputs the results to STDOUT.



### 2) Check overlap between (sub)samples ###

This script checks the overlap between one reference orthogroups set (the smallest one) and different subsamples of one large dataset to quantify the amount of orthogroups that are splitted when analysing larger datasets.

```
python overlap_checker.py
--ref <REFERENCE_ORTHOGROUP>
--ogs <ORTHOGROUP_FILES>
--out <OUTPUT_FOLDER>
```          

`--ref` specifies the OrthoFinder2 result while will serve as the reference orthogroup set.

`--ogs` specifies the OrthoFinder2 results which will be analyzed. This can be a single file or a comma-separated list of filenames.

`--out` specifies the result outputfile.



### 3) Check diversity of orthogroups ###

This script calculates Shannon diversity index and Eveness. A plot is generate to check for correlation between these values and the orthogroup size.


```
python diversity_check.py
--in <INPUT_OG_FILE>
--out <OUTPUT_FOLDER>
```          

`--in` specifies the OrthoFinder2 result will be analyzed.

`--out` specifies the output folder where all results will be stored.

							
							
### 4) Annotation of orthogroups ###

These scripts handle two important steps of the annotation process. In addition, a BLASTp search against the A. thaliana peptide sequence set is necessary.

```
python annotate_orthogroups1.py
--in <INPUT_FASTA_FILE_FOLDER>
--out <OUTPUT_FILE>
--og <ORTHOGROUP_FILE>
```          

`--in` specifies the folder with all underlying FASTA files. Sequences from all files will be loaded to construct a combined query FASTA file.

`--out` specifies the combined query FASTA file. Sequences will be stored in this file with the orthogroup name included as prefix in the sequence name.

`--og` specifies the OrthoFinder2 result file which contains the sequence IDs grouped into orthogroups.

							

#### running BLASTp search on compute cluster ###

```
python annotate_orthogroups2.py
--in <BLAST_RESULT_FILE>
--out <OUTPUT_FILE>
--og <ORTHOGROUP_FILE>
```          

`--in` specifies the BLAST result file. This can be reduced to the best hit per query to save computational resources.

`--out` specifies the result output file. The three best hits per orthogroup are stored.

`--og` specifies the OrthoFinder2 result file which contains the sequence IDs grouped into orthogroups.



### 5) Comparison of groups and identification of enriched orthogroups ###

Two contrasting groups can be compared to identify orthogroups that show overrepresentation of one of the groups. This can reveal lineage specific expansion and shrinking of gene families.

```
python enrichment_check.py
--in <BLAST_RESULT_FILE>
--out <OUTPUT_FOLDER>
--taxon <TAXON_FILE>

optional:
--anno <ORTHOGROUP_FILE>
```          

`--in` specifies an OrthoFinder2 result file which contains the sequence IDs in orthogroups.

`--out` specifies the output folder in which all result files will be saved. This folder will be created if it does not exist already.

`--taxon` specifies an the group membership of each taxon. This can be a simple table with taxon name in one column and group information in another column. Default: taxon name in second column and group information in fourth column.

`--anno` specifies an annotation file. The first column should contain the orthogroup ID and all following columns can contain annotation information. A TAB-separated file is expected. After extracting the ID, all columns starting with the second column will be merged into one annotation string. 





### References ###


