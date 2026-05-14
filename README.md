## Overview

This repository contains the computational pipeline and custom scripts developed for the genomic analysis of blood-derived Whole Genome Sequencing (WGS) data, as described in *Xeroderma pigmentosum monozygotic twins: same mutation divergent clinical outcomes*. 

The workflow is optimized for High-Performance Computing (HPC) environments using the SLURM workload manager and Singularity containers, ensuring reproducibility and scalability. It encompasses three main analytical modules:

1. **Somatic Variant Calling:** Automated SLURM job arrays utilizing GATK Mutect2 to accurately identify somatic mutations against a matched normal or Panel of Normals (PoN).
2. **Longitudinal Variant Filtering:** Custom Python and Bash scripts designed to isolate and filter temporally exclusive variants for each twin. This step could applies strict coverage (DP) and allele frequency (AF) thresholds to ensure high-confidence variant sets.
3. **Mutational Signatures Analysis:** Downstream analytical scripts to extract and characterize mutational profiles (e.g., SBS96 contexts), attributing the observed genomic patterns to specific underlying biological processes or deficiencies.
