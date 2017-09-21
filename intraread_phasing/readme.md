# intraread_phasing.py
Python script used to phase pairs of variants using individual reads that span both sites.

### Dependencies
- python
- pysam (https://github.com/pysam-developers/pysam)
- samtools (http://www.htslib.org/)


### Input

1. One or more **sam or bam files**

2. A **fasta file** for the reference genome

3. A **tab-delimited text file** with information about the variant pairs to be phased
    - This file should contain one line for each pair and each line should contain 6 columns in this order:
        1. Reference position for variant #1
        2. Reference allele for variant #1
        3. Alternative allele for variant #1
        4. Reference position for variant #2
        5. Reference allele for variant #2
        6. Alternative allele for variant #2

### Usage

To get get usage info:
```
intraread_phasing_v#.#.py  -h
```

Utilizing only reads with mapping quality ≥ 30 (-m) and where base quality at both variant sites ≥ 30 (-q):
```
intraread_phasing_v#.#.py -i tab-delim_pair_info.txt -o output.txt -q 30 -m 30 alignment.bam >stdout
```

### Output

    1. output.txt
        - Tab-delimited file with a header line ("File	Position 1	Position 2	Ref\Ref	Ref\Alt	Alt\Ref	Alt\Alt") and 1 row per bam/sam per variant pair
        - Each row will contain 7 columns:
            1. File: The sam/bam file analyzed
            2. Position 1: The reference position of the 1st variant
            3. Position 2: The reference position of the 2nd variant
            4. Ref\Ref: The # of reads with reference alleles at both positions
            5. Ref\Alt: The # of reads with the reference allele at position #1 and the alternative allele position #2
            6. Alt\Ref: The # of reads with the alternative allele at position #1 and the reference allele position #2
            7. Alt\Alt: The # of reads with alternative alleles at both positions
        -Only reads fitting into these four categories are reported in this table
        
    2. stdout
        - Raw counts of genotypes for each file, including genotypes that do not match the provided reference and alternative alleles
        - Use these counts to make sure that the analysis is working properly

### Options

  ```
Options:
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Text file with info for the position pairs to be
                        examined. One line per pair, tab-delimited with the
                        following columns: pos1, ref1, alt1, pos2, ref2, alt2
                        [None, OPT]
  -o OUT, --out=OUT     Name for output file, required if using the -i option.
                        [None, OPT/REQ]
  -q MINQUAL, --minQual=MINQUAL
                        Minimum base quals to be used in a geno [20]
  -m MAPQUAL, --mapQual=MAPQUAL
                        Minimum mapping quality for a read to be used [20]
  -1 FIRST, --first=FIRST
                        position in the reference of the left-most base to
                        phase
  -2 SECOND, --second=SECOND
                        position in the reference of the right-most base to
                        phase
  --R1only              Use this flag to only consider R1
  --R2only              Use this flag to only consider R2
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
