# plot_stranded_coverage.py
Python script used to generate strand-specific depth of coverage plots


### Dependencies
- python
- pysam (https://github.com/pysam-developers/pysam)
- samtools (http://www.htslib.org/)
- matplotlib (https://matplotlib.org)

### Input

1. One or more **sam or bam files** specified as arguments on the command line

### Usage

To get get usage info:
```
plot_stranded_coverage_v#.#.py  -h
```

Generate individual coverage plots for three different samples. Utilizing only R1 reads with mapping quality ≥ 30 (-q 30):
```
plot_stranded_coverage_v#.#.py -q 30 sample1.bam sample2.bam sample3.bam
```

Generate individual coverage plots for three different samples. Utilizing **only R2 reads (--useR2)** with mapping quality ≥ 30 (-q 30):
```
plot_stranded_coverage_v#.#.py -q 30 --useR2 sample1.bam sample2.bam sample3.bam
```

Generate individual coverage plots for three different samples. Utilizing only R1reads with mapping quality ≥ 30 (-q 30):

**Smooth coverage plots using a 500 nucleotide sliding window, with a step size of 50 nucleotides (-s 500,50)**
```
plot_stranded_coverage_v#.#.py -q 30 -s 500,50 sample1.bam sample2.bam sample3.bam
```


### Output

    Note: 
    - If using R1, Forward reads are from RNA molecules that were the reverse complement of the reference (E.g., if using positive strand reference, R1 Forward reads = negative strand)
    - If using R2, Forward reads are from RNA molecules that were the same strand as the reference (E.g., if using positive strand reference, R2 Forward reads = positive strand)

    1. sample1.bam_ref_strandcov.pdf, sample2.bam_ref_strandcov.pdf, sample3.bam_ref_strandcov.pdf
        - pdf coverage plots for each individual sample
        - Each plot will contain two lines:
            - Red = Forward strand coveage
            - Green = Reverse strand coverage
        - ref = name of reference sequence used for alignment

### Options

  ```
Options:
  -h, --help            show this help message and exit
  --useR2               Use this flag if you want use R2 instead of R1 [false]
  -q MAPQ, --mapq=MAPQ  minimum mapping quality to be used [20]
  -s SMOOTH, --smooth=SMOOTH
                        Use this option if you want to smooth the coverage
                        plots. Specify window size and slide, comma separated
                        [None]
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
