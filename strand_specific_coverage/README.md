# strand_specific_coverage

A collection of python scripts used for analyzing stranded RNAseq data. Designed for use with small viral genomes. 

    1. strand_ratio_counts.py
        - Starting with a sam/bam file, calulates the # and proportion of reads mapping to each genomic strand
    2. plot_stranded_coverage.py
        - Starting with a sam/bam file, generates a coverage plot with separate lines for each genomic strand.
    3. plot_stranded_coverage_multinorm.py
        - Starting with multiple sam/bam files aligned to the same reference, generates normalized, strand-specific coverage plots
        - Normalized coverage calculations are focused on user-specified windows (e.g., ORFs)
    4. strandspec_covplot.py
        - Starting with multiple sam/bam files aligned to the same reference, generates normalized and non-normalized, strand-specific coverage plots
        - Coverage plots generated using a sliding window across the entire genome


Copyright (C) 2017  Jason Ladner

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with each program.  If not, see <http://www.gnu.org/licenses/>.
