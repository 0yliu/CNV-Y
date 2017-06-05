# denovo CNV detection
De novo Duplication source detection

Before starting processing the data, there are several cleanning steps needed to be done.
1. Check the family structures, trios, quads, etc. 
2. Check the tissue types for each samples. For example, in Craniofacial data set, there are blood samples (mainly), buccal samples, saliva samples, and samples being WGA (whole genome amplification) because lacking of enough DNA. In Autism data set, a certain proportion of samples are cell line samples, which we need to be cautious about.
3. Check the genotype data quality (mainly the BAF and the LRR, almost 1000 samples in AGP dataset has no information about BAF and LRR, we remove them from it), as well as the call rate. If the missing call rate of a sample is larger than 5%, then it might be removed.

