# ID_revert_TC_clusters.py
Python script used to identify and revert clusters of T-to-C mutations indicative of ADAR hyperediting

### Dependencies
- python ≥2.6


### Input

1. An **aligned fasta file**
2. **Name of sequence** in the alignment to be used as a reference. Ideally this will be an ancestral sequence, will be used to infer directionality of substitutions.
3. **Size of sliding window** to be used for identifying clusters
4. **Minimum number of T-to-C changes** within a window to be flagged as a cluster and for the substitutions to be reverted to Ts

### Usage

To get get usage info:
```
ID_revert_TC_clusters.py  -h
```

Identify and revert any clusters of at least 4 T-to-C changes within a 200 nucleotide sliding window:
```
ID_revert_TC_clusters.py -f aligned.fasta -r refname -w 200 -n 4
```
Note: By default, the name of the fasta file will be used as the basename for the output files. If you want a different base name, use "-o"


### Output

    1. aligned_ADARrversed.fasta
        - Your aligned fasta file, but with identified clusters of T-to-C changes reverted to Ts

    2. aligned_clusters.txt
        - Tab-delimited text file with info about the identified clusters
        - Overlapping clusters are merged for reporting purposes
        - Column 1 = name of sequence from fasta with T-to-C cluster
        - Column 2 = comma-delimited list of nucleotide positions within the cluster

    3. aligned_refname.vcf
        - vcf format file summarizing substitutions that differentiate each sequence from the chosen reference
    
    4. Information about the sequence context (upstream and downstream nucleotides) of all putative ADAR-related changes are printed to standard error

### Options
```
Options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta=FASTA
                        Aligned fasta. [None, REQ]
  -r REF, --ref=REF     Name of sequence in the alignment to use as a
                        reference. [None, REQ]
  -o OUT, --out=OUT     Base name for output files. [None, OPT]
  -g GAP, --gap=GAP     Character used for gaps. [-]
  -v VCF, --vcf=VCF     vcf file. [None, OPT]
  -w WIN, --win=WIN     Window size to check for enrichment. [200]
  -n NUM, --num=NUM     Min number of changes for enrichment. [4]
```


Copyright (C) 2017  Jason Ladner

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
