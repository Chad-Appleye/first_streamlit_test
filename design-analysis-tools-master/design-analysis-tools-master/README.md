# design-analysis-tools
Tools for analysis of genomes and CoPrimers in assay development.

## primer_mismatch.py
A module for analyzing the mismatches for a given set of primers

### Usage
Example for a forward and reverse CoPrimer design. Also accepts normal primers and probes.

    from primer_mismatch import parse_alignment, PrimerCoverage

    aln = parse_alignment('path_to_alignment_file')
    forward_primers = ['tGGAGGATGAAGAAGATG[BHQ-2][Spacer 18]aaaaaaaaa[Spacer 18]tgATGATCTTACAG']
    reverse_primers = ['AGTTTCAGCTGCTCG[Spacer 18]aaaaaaaaa[Spacer 18]gaCTCCCACCG']
    forward_gaps = [3]  # gap is necessary to find position in alignment
    reverse_gaps = [4]
    coverage = PrimerCoverage(aln, forward_primers, reverse_primers, 
        forward_gaps=forward_gaps, reverse_gaps=reverse_gaps)
    coverage.display_report()
    
### Pip install
`pip install -e git+https://github.com/Co-Diagnostics/design_analysis_tools.git@master#egg=design-analysis-tools`

