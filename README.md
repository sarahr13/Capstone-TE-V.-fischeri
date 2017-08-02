# Capstone-TE-V.-fischeri
Code and data for STAT 5010-5020 project (Fall '16 -- Spring '17)

Client: Eric Stabb, Department of Microbiology, University of Georgia.

#### Introduction
Transposable elements (TEs) are pieces of DNA that can "jump" around in the genome. TEs can be used as a mutagenesis tool if it inserts in the desired location. Our client, Eric Stabb, is interested in whether they have an insertion bias in the bacteria, Vibrio fischeri. It is believed that the insertions are not random. Our data consist of TE insertions in 9542 V. fischeri mutants. The TE is noted to be on chromosome 1, chromosome 2, or the plasmid. We investigate the insertion preference using AT/GC-content and the coding vs. non-coding percentage of the insertions.

#### Methods
- Partition the genome into "bins" and calculate the number of insertions, GC-content and coding percentage of each bin
  - Chromosome bin sizes: 1000 bp, 5000 bp, 10000 bp
  - Plasmid bin sizes: 750 bp, 1000 bp, 1500 bp
- Test for uniformity of insertions *across* the entire genome (preference for either chromosome or plasmid) 
- Test for uniformity of insertions *within* each genetic element (specific locations within each genetic element)
- Simulate insertion data under the hypothesis of uniformly distributed insertion locations 

#### Results
###### Uniformity tests
- Rejected null hypothesis of uniform distribution of insertion locations for both uniformity tests

###### Negative binomial regression
For chromosome 1 with a bin width of 5000 bp, we found that GC percent, coding percent, and the interaction between these two characteristicswere all highly significant predictors.  For chromosome 2, only GC content was significant.  None of the predictors were significant for the plasmid at a bin width of 1000 bp.

#### Conclusion
Based on our bin widths, we found that GC-content, coding percentage, and their interactionwere all significant in chromosome 1. For chromosome 2, only GC-content was not significant, exceptfor when we set the bin width to 1000 bp.  We did not find GC-content, coding percentage, or theirinteraction to be significant for the plasmid at any tested bin width.  This is likely due to the smallsample size of the plasmid.  We also identified specific regions within the each genetic element thatcontain significantly more insertions than expected under a random insertion distribution
