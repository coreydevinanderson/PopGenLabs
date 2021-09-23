This lab was developed to teach undergraduate students basic population genetic calculations and tests, based on analysis of real microsatellite DNA data. This was designed as a two week lab (i.e. two 3-hour sessions), with homework each week. The students use the armadillo data in class, and then repeat the analyses for a different data set (grasshopper sparrows).

There is a strong MS Excel componenet to this lab and, after doing calculations manually (in Excel), they use software to facilitate calculations and do basic tests. Historically, we used the Excel Microsatellite Toolkit to do basic summary statistics and format data for other software. We conducted tests of significance of the system of mating inbreeding coefficient in FSTAT and tests of genotypic linkage disequilibrium in GENEPOP. However, the Excel Microsatellite Toolkit is no longer available (albeit it still works, with some effort) and FSTAT is defunct. Hence, I am in the process of replacing these external tools with R code. So far, I have code to do basic summaries of allele and genotype frequencies, and to calculate and conduct permutation tests of the system of mating inbreeding coefficient. The original lab instructions have not yet been updated to include the new R code.

There are four R files that contain different functions:

* read_files.R
- This file contains basic R code for reading in the comma delimited files (from the working directory) and an additional function (summary_df()) that provides a basic summary of the imported file. The microsatellite data must be in diploid two column format, with the first column being IDs.

* allele_summary.R
- The output 'geno' data.frame from read_files.R is used by the allele_summary() function. This function outputs a basic summary of alleles for a specifed locus.

* genotype_summary.R
- The output 'geno' data.frame from read_files.R is used by the genotype_summary() function. This function outputs a basic summary of the genotypes for a specified locus.

* het_and_f.R
- This file contains six functions (f(), f_multi(), f_avg(), f_test(), f_test_multi(), and f_test_avg()) for calculating and testing the system of mating inbreeding coefficient. These functions have dependencies on the output of genotypes_summary.R and allele_summary.R, as well as on eachother.

* R_easy_peasy.R
- This file contains all the functions. If you run the entire block of code, it will output two comma delimited files into the working directory with the results from f_test_multi() and f_test_avg(). This is a nice option if you don't want to teach the students about R and just want a result.

The resources file contains the armadillo and genotype data, the original lab instructions, and an additional .xlsx file with the uDNA data on separate sheets.

Please email me if you would like the solutions.

coreanderson@valdosta.edu
  
