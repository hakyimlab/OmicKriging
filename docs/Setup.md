## Depedencies

Ensure the following R packages are installed on your environment. OmicKriging utilises these packages.

- ROCR
- doParallel
- parallel
- foreach
- iterators

**NB:** These libraries are required by the OmicKriging.

## Phenotype file

The phenotype file used in your analysis should;
- Contain phenotypes on the columns and samples on the rows
- Not contain NA values, these break the analysis pipeline.
   - Mutate the NA values with column averages for each phenotype.
- The individual IDs used in the phenotype file should be same as those in the grm file.
- Contain an equal number or a subset of the individuals in the grm matrix.
  - `E.g` if a grm file has 997 individuals the phenotype should have 997 individuals or less.
  - If the phenotype file contains many individuals than those in the grm it will break the analysis pipeline.
  
- The sample ids column should be denoted as `IID`.
- The phenotype file should contain the sample id `IID` and family id `FID` in the first two columns. Phenotypes should start from column 3 to the end.
