Pipelines designed to run two-sample medelian randomization and coloc between GWAS summary statistics and molecular QTL summary statistics generated from fastQTL/tensorQTL.

The coloc pipeline requires:
- The GWAS summary statistics
- A table of GWAS significant independent signals
- The full molecular QTL summary statisticsfrom in fastQTL/tensoQTL format
- The result from the permutation analysis from fastQTL/tensoQTL
- A reference genotype file

For Coloc, the following steps are performed:
- the GWAS summary statistics are converted to a standardised, tabix indexed, VCF file format. This step includes liftover between genome builds.
- The full summary statistics of the molecular QTLs are filtered to only keep those where the lead SNP is +/- 1mb from the GWAS lead SNP.
- The GWAS and molecular QTL summay statistics are harmonised according to the reference genotype. This can be the in-sample genotype from the GWAS or molQTL or an external one but, in order to run coloc.SuSiE, it needs to be the same that will be used to generate the LD matrix.
- Colocalization with coloc.abf for the overlapping regions.
- For regions with little evidence of colocalization, coloc.SuSiE is used to test for colocalization on multiple causal variants.

TODO: Always take care of confilcts ! We will hear about that now !
Project Lead: Norbert
