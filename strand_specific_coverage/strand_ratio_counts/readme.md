# strand_ratio_counts.py
Python script used to calculate the proportion of reads originating from +/- sense RNA molecules

### Dependencies
- python
- pysam (https://github.com/pysam-developers/pysam)
- samtools (http://www.htslib.org/)


### Input

3. A **tab-delimited text file** specifying the **bam/sam files** to be analyzed 
    - This file should contain a header line (header values are arbitrary) and at least 3 columns (only the 3rd is utilized):
        1. Sample info of user's choice (e.g., sample name)
        2. Sample info of user's choice (e.g., sample type)
        3. Name of bam/sam file (can be absolute or relative path)

### Usage

To get get usage info:
```
strand_ratio_counts_v#.#.py  -h
```

Default usage, considering reads mapped anywhere in the reference:
```
strand_ratio_counts_v#.#.py -i tab-delim_bam_info.txt -o output.txt
```

Only considering reads that overlap by at least 50 nucleotides (-v 50) with the first 2000 nucleotides of the reference (-b 0 -e 1999):
```
strand_ratio_counts_v#.#.py -i tab-delim_bam_info.txt -o output.txt -v 50 -b 0 -e 1999
```


### Output

    Notes: 
    - Only counts R1 and only if part of 'proper pair'
    - Forward R1 reads are from RNA molecules that were the reverse complement of the reference (E.g., if using positive strand reference, R1 Forward reads = negative strand)


    1. output.txt
        - Tab-delimited file with a header line per bam/sam (identical to input file, but with 5 additional columns)
        - Each row will contain 8 columns:
            1. 1st column from the input file
            2. 2nd column from the input file
            3. 3rd column from the input file (sam/bam names)
            4. Proportion of forward reads (-99 if no reads mapping)
            5. Proportion of reverse reads (-99 if no reads mapping)
            6. Strand ratio: #reverse/#forward (-99 if no forward reads)
            7. Number of forward reads
            8. Number of reverse reads

### Options

  ```
Options:
  -h, --help            show this help message and exit
  -i INP, --inp=INP     Name for input file [None]
  -o OUT, --out=OUT     Name for output file [None]
  -b BEG, --beg=BEG     1st base in range
  -e END, --end=END     last base in range
  -v OVERLAP, --overlap=OVERLAP
                        required overlap for reads to be inlcuded [10]
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
