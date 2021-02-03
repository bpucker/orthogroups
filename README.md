# orthogroups
script related to the filtering and analysis of orthogroups in the Caryophyllales


### 1) Remove orthogropus with low taxa represenation (likely contamination/artifacts) ###

This script takes the orthogroup.tsv file produced by OrthoFinder2 and transfers only orthogroup with a certain level of taxa representation into a new file. The --optimize options evaluates different cutoff values and supports the identification of a suitable one.

filter_OG_by_taxa.py
							--out <OUTPUT_OG_FILE>
							--in <INPUT_OG_FILE>
							--taxa <MINIMAL_TAXA_NUMBER> | --optimize


### 2) Check overlap between (sub)samples ###

This script checks the overlap between one reference orthogroups set (the smallest one) and different subsamples of one large dataset to quantify the amount of orthogroups that are splitted when analysing larger datasets.

python overlap_checker.py
							--ref <REFERENCE_ORTHOGROPU>
							--ogs <ORTHOGROUP_FILES>
							--out <OUTPUT_FOLDER>



### References ###


