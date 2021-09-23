## genotype_summary() ##

# Dependencies:  must first  create the data.frame (= geno)

# Description: Summarizes genotype data for a locus

## Arguments:

# your_df: a data.frame with your genotype data (= geno)
# col_1, col_2: the columns with your uDNA data (in diploid two-column format)

# Returns: the count (freq) for each genotype and the relative freq


## Run the function first (from # Start to # End), then call it.

# Start
genotype_summary <- function(your_df, col_1, col_2){

  col_x <- your_df[, col_1]
  col_y <- your_df[, col_2]

  col_x[col_x == 0] <- NA
  col_y[col_y == 0] <- NA
	
  geno_df <- data.frame(col_x, col_y)
  geno_df <- geno_df[order(col_x, col_y),]
  missing <- sum(is.na(geno_df$col_x))  # number of missing genotypes
  geno_df <- na.omit(geno_df)
  geno_unique <- unique(geno_df)

  geno_unique_comb <- paste(geno_unique[, 1], geno_unique[, 2], sep = "")
  geno_df_comb <- paste(geno_df[, 1], geno_df[, 2], sep = "")

  counts <- integer(dim(geno_unique)[1]) # make empty vector
  for (i in 1:length(geno_unique_comb)){
    counts[i] <- sum(geno_unique_comb[i] == geno_df_comb)
  }

  geno_unique["freq"] <- counts
  geno_unique["r_freq"] <- round(counts / sum(counts), 4)
  row.names(geno_unique) <- c()
	
return(geno_unique)
}
# End

## Change the column numbers as required

genotype_summary(geno, 2, 3)

genotype_table <- genotype_summary(geno, 2, 3)  # change col numbers as required

print(genotype_table)
