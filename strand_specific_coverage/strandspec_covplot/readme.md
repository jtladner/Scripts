# strandspec_covplot.py
Python script used to strand-specific coverage plots.
    - Quick and dirty plots are generated with matplotlib
    - However, the raw data for the plots is also output and can be migrated to other graphing softwrae options

### Dependencies
- python
- pysam (https://github.com/pysam-developers/pysam)
- samtools (http://www.htslib.org/)
- numpy (www.numpy.org)
- matplotlib (https://matplotlib.org)

### Input

1. One or more **sam or bam files** specified as arguments on the command line
    - The reference needs to be the same in all specified alignment files
    
### Usage

To get get usage info:
```
strandspec_covplot_v#.#.py  -h
```

Analyzing 3 different bams, plotting coverage using only R1 (default) and utilizing only reads with mapping quality ≥ 30 (-q):
```
strandspec_covplot_v#.#.py -o out -q 30 sample1.bam sample2.bam sample3.bam
```

Analyzing 3 different bams, plotting coverage using only R1 (default) and utilizing only reads with mapping quality ≥ 30 (-q):

**Coverage plots smoothed using a sliding window of 500 nt and a slide of 100 nt (--smooth 500,100)**
```
strandspec_covplot_v#.#.py -o out -q 30 --smooth 500,100 sample1.bam sample2.bam sample3.bam
```

Analyzing 3 different bams, plotting coverage using only **R2 (--useR2)** and utilizing only reads with mapping quality ≥ 30 (-q):
```
strandspec_covplot_v#.#.py -o out -q 30 --useR2 sample1.bam sample2.bam sample3.bam
```

Analyzing 3 different bams, plotting coverage using R1 (default) and **unpaired reads (--useUnpaired)**, Utilizing only reads with mapping quality ≥ 30 (-q):
```
strandspec_covplot_v#.#.py -o out -q 30 --useUnpaired sample1.bam sample2.bam sample3.bam
```

Analyzing 3 different bams, plotting coverage using only R1 (default) and utilizing only reads with mapping quality ≥ 30 (-q):

**Only generating combined coverage plots (--noIndivPlots)**
```
strandspec_covplot_v#.#.py -o out -q 30 --noIndivPlots sample1.bam sample2.bam sample3.bam
```


### Output

    Note: 
    - If using R1, Forward reads are from RNA molecules that were the reverse complement of the reference (E.g., if using positive strand reference, R1 Forward reads = negative strand)
    - If using R2, Forward reads are from RNA molecules that were the same strand as the reference (E.g., if using positive strand reference, R2 Forward reads = positive strand)


    1. out_rawcov.txt
        - Tab-delimited text file (5 columns w/ header) containing the raw coverage information used to generate the plots
            1. sam/bam file name
            2. Reference sequence name
            3. Nucleotid position (1-based)
            4. Forward strand coverage
            5. Reverse strand coverage
            
    2. out_normcov.txt
        - Tab-delimited text file  (5 columns w/ header) containing the normalized coverage information used to generate the plots
            1. sam/bam file name
            2. Reference sequence name
            3. Nucleotid position (1-based)
            4. Normalized forward strand coverage
            5. Normalized reverse strand coverage

    3. out_ref_For_combo.pdf, out_ref_Rev_combo.pdf
        - Separate combo plots for forward and reverse strand reads
            - Solid line = normalized mean across samples
            - Ribbon = standard deviation across samples
        - ref = name of reference sequence used for alignment

    4. sample1.bam_ref_rawcov.pdf,sample1.bam_ref_normcov.pdf, ...
        - Coverage plots, both for raw coverage and normalized coverage, for each individual sample
        - Each plot will contain two lines:
            - Red = Forward strand coveage
            - Green = Reverse strand coverage
        - Will **NOT** be produced if --noIndivPlots flag is used
        - ref = name of reference sequence used for alignment

### Options

  ```
Options:
  -h, --help            show this help message and exit
  -o OUT, --out=OUT     Base name for output files [None]
  -s SMOOTH, --smooth=SMOOTH
                        Use this option if you want to smooth the coverage
                        plots. Specify window size and slide, comma separated
                        [None]
  --useR2               Use this flag if you want use R2 instead of R1 [false]
  --useUnpaired         Use this flag if you want use reads not specified as
                        R1 or R2 [false]
  -q MAPQ, --mapq=MAPQ  minimum mapping quality to be used [20]
  --noIndivPlots        Use this flag if you do not want to create plots for
                        individual bams [false]
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
