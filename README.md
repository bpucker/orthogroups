# orthogroups
script related to the filtering and analysis of orthogroups in the Caryophyllales


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

python overlap_checker.py
							--ref <REFERENCE_ORTHOGROPU>
							--ogs <ORTHOGROUP_FILES>
							--out <OUTPUT_FOLDER>




### 3) Check diversity of orthogroups ###

This script calculates Shannon diversity index and Eveness. A plot is generate to check for correlation between these values and the orthogroup size.

python diversity_check.py
							--out <OUTPUT_FOLDER>
							--in <INPUT_OG_FILE>
							
							
### 4) Annotation of orthogroups ###

These scripts handle two important steps of the annotation process. In addition, a BLASTp search against the A. thaliana peptide sequence set is necessary.

python annotate_orthogroups1.py
							--out <OUTPUT_FILE>
							--og <ORTHOGROUP_FILE>
							--in <INPUT_FASTA_FILE_FOLDER>
							

RUN BLASTp search on compute cluster

python annotate_orthogroups2.py
							--out <OUTPUT_FILE>
							--og <ORTHOGROUP_FILE>
							--in <BLAST_RESULT_FILE>
							feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de





### References ###


