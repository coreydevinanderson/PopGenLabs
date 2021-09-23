
### allele_summary() ###

# Dependencies: must first create the data.frame (= geno)

# Description: Summarizes allele data for a specified locus.

## Arguments:

# your_df: a data.frame with your genotype data (= geno); see read_files.R
# col_1, col_2: the columns associated with the locus you want to summarize...
# ...in diploid two-column format.

# Returns: N2, N, count of missing data, alleles, allele freqs, & relative freqs

##########

# Run the function first (from # Start to # End), then call it.

# Start
allele_summary <- function(your_df, col_1, col_2){

  col_x <- your_df[, col_1]
  col_y <- your_df[, col_2]

  col_x[col_x == 0] <- NA
  col_y[col_y == 0] <- NA
  
  N_raw <- length(c(col_x, col_y))

  all_alleles <- as.factor(c(as.character(na.omit(col_x)), as.character(na.omit(col_y))))

  N_2 <- length(all_alleles)
  N <- N_2 / 2
  unique_alleles <- levels(all_alleles)

  allele_freqs <- table(all_alleles)
  allele_r_freqs <- round(table(all_alleles)/sum(allele_freqs), 4)

return(list(N2 = sum(allele_freqs), 
		N = sum(allele_freqs) / 2,
		missing = (N_raw / 2) -(sum(allele_freqs) / 2), 
		alleles = unique_alleles, 
		allele_freqs = allele_freqs, 
		allele_r_freqs = allele_r_freqs))
}
# End

## Change the column numbers as required

allele_summary(geno, 2, 3)

allele_list <- allele_summary(geno, 2, 3)  # this saves the output as an allele_table

print(allele_list)  # if you just want to view it





	



