### Functions to calculate the system of mating inbreeding coefficient for a single population.
### Based on observed and expected heterozygosity.

## Functions:

## f()  # calculate Hobs, Hexp, and f for a single locus
## f_multi()  # Summarize heterozygosity and f for all loci in your data.frame
## f_avg()  # Average f over all loci in sample
## f_test() # Perform a permutation test of the significance of f for a single locus.
## f_test_multi() # Perform a permutation test of the significance of f for all loci.
## f_test_avg() # Performa a permutation test of the average f over all loci.

########################################################################

## f() ##

# Description: Calculates Hobs, Hexp, and f for a single locus

# Dependencies: Requires the output list from allele_summary() and the output 
# ...data.frame from genotype_summary()

# For example: 
# genotype_table <- allele_summary(geno, 2, 3)
# allele_list <- genotype_summary(geno, 2, 3)
# Now you can call the function:
# f(genotype_table, allele_list)

## Arguments:

# genotype_data: data.frame obejct from genotype_summary()
# allele_data: list object from allele_summary()

# Returns Hobs, Hexp, Hexp_nei, f_raw, f_corrected 


##Run the function first (from # Start to # End), then call it.

# Start
f <- function(genotypes_out, alleles_out){
	
  N <- sum(genotypes_out[, 3])	
  t_f <- genotypes_out[, 1] == genotypes_out[, 2]
  homo_obs <- sum(t_f * genotypes_out[, 3]) / N
  het_obs <- 1 - homo_obs

  allele_r_freqs <- unname(alleles_out$allele_r_freqs)
  het_exp <- 1 - sum(allele_r_freqs ^ 2)
  N2 <- alleles_out$N2
  Hexp_nei <- (N2 / (N2 - 1)) * het_exp
  
return(list(Hobs = het_obs,
		Hexp = het_exp,
	      Hexp_nei = Hexp_nei,
		f_raw = round(1 - (het_obs / het_exp), 3),
		f_corrected = round(1 - (het_obs / Hexp_nei), 3)))
}
# End

f(genotype_table, allele_list)


########################################################################

## f_multi ##

# Description: Summarizes heterozygosity and f for all loci in your data.frame

# Dependencies: must first create the data.frame (= geno) and that you run...
# ...the code for f() prior to use.

## Arguments:

# your_df: a data.frame with your genotype data (= geno)

# Returns: a data.frame with Hobs, Hexp, Hexp_nei, f_raw, f_corrected for...
# ...each locus. 


## Run the function first (from # Start to # End), then call it.

# Start
f_multi <- function(your_df){
  
  Hobs <- c()
  Hexp <- c()
  Hexp_nei <- c()
  f_raw <- c()
  f_corrected <- c()
	
  i <- 2

  limit <- dim(your_df)[2]+1
	
  while(i < limit){
    alleles_temp <- allele_summary(your_df, i, i+1)
    genotypes_temp <- genotype_summary(your_df, i, i+1)
    info <- f(genotypes_temp, alleles_temp)
    Hobs <- append(Hobs, info$Hobs)
    Hexp <- append(Hexp, info$Hexp)
    Hexp_nei <- append(Hexp_nei, info$Hexp_nei)
    f_raw <- append(f_raw, info$f_raw)
    f_corrected <- append(f_corrected, info$f_corrected) 
    i <- i + 2
  }

return(data.frame(Hobs, Hexp, Hexp_nei, f_raw, f_corrected))
}
# End

f_multi(geno)

f_table <- f_multi(geno)  #if you want to store the output object for f_avg

write.csv(f_table, "f_table.csv")  # If you want to output a table to Excel
                                   # Output table will be in your working directory  


########################################################################


## f_avg ##

# Description: averages the system of mating inbreeding coefficient over all...
# ...loci in your sample.

# Dependencies: Requires output data.frame from f_multi()

## Arguments:

# fs: a data.frame from f_multi() containing raw and corrected f values


## Run the function first (from # Start to # End), then call it.

# Start
f_avg <- function(fs){
	mean_f_raw <- mean(fs$f_raw)
	mean_f_corr <- mean(fs$f_corrected)
	
return(c("mean f_raw"= mean_f_raw, "mean f_corr" = mean_f_corr))
}
# End

f_avg(f_table)

write.csv(f_avg, "f_avg.csv")   # If you want to output a table to Excel
                                # Output table will be in your working directory


#########################################################################


## f_test ##

# Description: conducts a permutation test for f.

# Dependencies: must first create the data.frame (= geno) and must run code...
# ...for allele_summary() and genotype_summary()

## Arguments:

# your_df: a data.frame with your genotype data (= geno)
# col_1, col_2: the columns with your uDNA data (in diploid two-column format)
# nsims: the number of permutations

# Returns: a named vector containing the estimate of f and the P_value (one-tailed)


## Run the function first (from # Start to # End), then call it.

# Start
f_test <- function(your_df, col_1, col_2, nsims){
  
  genotype_summary2 <- function(your_df, col_1, col_2){

    col_x <- your_df[, col_1]
    col_y <- your_df[, col_2]

    col_x[col_x == 0] <- NA
    col_y[col_y == 0] <- NA
	
    geno_df <- data.frame(col_x, col_y)
    geno_df <- geno_df[order(col_x, col_y),]
    missing <- sum(is.na(geno_df$col_x))  # number of missing genotypes
    geno_df <- na.omit(geno_df)
  
    geno_dims <- dim(geno_df)
    allele_rvect <- sample(c(geno_df[, 1], geno_df[, 2]))
    col_xx <- allele_rvect[1:geno_dims[1]]  
    col_yy <- allele_rvect[(geno_dims[1] + 1):(geno_dims[1] * geno_dims[2])]
    geno_df2 <- data.frame(col_xx, col_yy)
    geno_df2 <- geno_df2[order(col_xx, col_yy),]  

    geno_unique <- unique(geno_df2)

    geno_unique_comb <- paste(geno_unique[, 1], geno_unique[, 2], sep = "")
    geno_df_comb <- paste(geno_df2[, 1], geno_df2[, 2], sep = "")

    counts <- integer(dim(geno_unique)[1]) # make empty vector
    for (i in 1:length(geno_unique_comb)){
      counts[i] <- sum(geno_unique_comb[i] == geno_df_comb)
    }

    geno_unique["freq"] <- counts
    geno_unique["r_freq"] <- round(counts / sum(counts), 4)
    row.names(geno_unique) <- c()
	
  return(geno_unique)
  }

  alleles_out <- allele_summary(your_df, col_1, col_2)
  genotypes_out <- genotype_summary(your_df, col_1, col_2)
	
  N <- sum(genotypes_out[, 3])	
  t_f <- genotypes_out[, 1] == genotypes_out[, 2]
  homo_obs <- sum(t_f * genotypes_out[, 3]) / N
  het_obs <- 1 - homo_obs

  allele_r_freqs <- unname(alleles_out$allele_r_freqs)
  het_exp <- 1 - sum(allele_r_freqs ^ 2)
  N2 <- alleles_out$N2
  Hexp_nei <- (N2 / (N2 - 1)) * het_exp
  f_obs = round(1 - (het_obs / Hexp_nei), 3)

  f_permuted <- numeric(nsims)

  for(i in 1:nsims){
    genotypes_rand <- genotype_summary2(your_df, col_1, col_2)	
    t_f_rand <- genotypes_rand[, 1] == genotypes_rand[, 2]
    homo_obs_rand <- sum(t_f_rand * genotypes_rand[, 3]) / N
    het_obs_rand <- 1 - homo_obs_rand
    f_permuted[i] <- round(1 - (het_obs_rand / Hexp_nei), 3)
  }
    
  f_permuted <- append(f_permuted, f_obs)

  if(f_obs >= 0){  
   P_val <- sum(f_permuted >= f_obs) / length(f_permuted)
  } else {
   P_val <- sum(f_permuted < f_obs) / length(f_permuted)
  }

  out_vect <- c(f_corrected = f_obs, P_value = P_val)

return(round(out_vect, 3))
}
# End

f_test(geno, 2, 3, nsims = 9999)


#########################################################################

## f_test_multi ##

# Description: Conducts test of f for each locus.

# Dependencies: must first create the data.frame (= geno) and  run code...
# ...for f_test() (and its required functions) prior to use.

## Arguments:

# your_df: a data.frame with your genotype data (= geno)
# nsims: the number of permutations

# Returns: a matrix containing f for each locus and P-values from the permuation tests.


## Run the function first (from # Start to # End), then call it.

# Start
f_test_multi <- function(your_df, nsims){

  i <- 2
  limit <- dim(your_df)[2]+1
  out_matrix <- c()
	
  while(i < limit){
    out_matrix <- rbind(out_matrix, f_test(your_df, i, i + 1, nsims)) 
    i = i + 2
  }


  # Use locus names as row names

  locus_names_raw <- names(geno)
  locus_trimmed <- c()
  limit2 <- dim(geno)[2]
  count <- 2

  repeat{
   locus_trimmed <- append(locus_trimmed, locus_names_raw[count])
   count <- count + 2
   if(count > limit2){
     break
   }
  }

  locus_names <- c()
 
  for(j in 1:length(locus_trimmed)){
    locus_names <- append(locus_names, substr(locus_trimmed[j], 1, nchar(locus_trimmed[j]) - 1))
  }
  
  rownames(out_matrix) <- locus_names
  
return(out_matrix)
}
# End

f_test_multi(geno, 9999)

f_tests <- f_test_multi(geno, 9999)

write.csv(f_tests, "f_tests.csv")  # If you want to output a table to Excel
                                   # Output table will be in your working directory


#########################################################################


## f_test_avg ##

# Description: conducts a permutation test for the mean of f over all loci

# Dependencies: must first create the data.frame (= geno) and  run code...
# ...for allele_summary(), genotype_summary(), and f_multi() (and its...
# ...associated function) prior to use.

## Arguments:

# your_df: a data.frame with your genotype data (= geno)
# nsims: the number of permutations

# Returns: a named vector containing f_avg and the P-value.


## Run the function first (from # Start to # End), then call it.

# Start
f_test_avg <- function(your_df, nsims){

  genotype_summary2 <- function(your_df, col_1, col_2){

    col_x <- your_df[, col_1]
    col_y <- your_df[, col_2]

    col_x[col_x == 0] <- NA
    col_y[col_y == 0] <- NA
	
    geno_df <- data.frame(col_x, col_y)
    geno_df <- geno_df[order(col_x, col_y),]
    missing <- sum(is.na(geno_df$col_x))  # number of missing genotypes
    geno_df <- na.omit(geno_df)
  
    geno_dims <- dim(geno_df)
    allele_rvect <- sample(c(geno_df[, 1], geno_df[, 2]))
    col_xx <- allele_rvect[1:geno_dims[1]]  
    col_yy <- allele_rvect[(geno_dims[1] + 1):(geno_dims[1] * geno_dims[2])]
    geno_df2 <- data.frame(col_xx, col_yy)
    geno_df2 <- geno_df2[order(col_xx, col_yy),]  

    geno_unique <- unique(geno_df2)

    geno_unique_comb <- paste(geno_unique[, 1], geno_unique[, 2], sep = "")
    geno_df_comb <- paste(geno_df2[, 1], geno_df2[, 2], sep = "")

    counts <- integer(dim(geno_unique)[1]) # make empty vector
    for (i in 1:length(geno_unique_comb)){
      counts[i] <- sum(geno_unique_comb[i] == geno_df_comb)
    }

    geno_unique["freq"] <- counts
    geno_unique["r_freq"] <- round(counts / sum(counts), 4)
    row.names(geno_unique) <- c()
	
  return(geno_unique)
  }

  f_avg_vect <- numeric(nsims) 

  for(j in 1:nsims){ 

    k <- 2
    limit <- dim(your_df)[2]+1
    f_vect <- c()
	
    while(k < limit){

      alleles_out <- allele_summary(your_df, k, k + 1)
      genotypes_out <- genotype_summary(your_df, k, k + 1)
	
      N <- sum(genotypes_out[, 3])	

      allele_r_freqs <- unname(alleles_out$allele_r_freqs)
      het_exp <- 1 - sum(allele_r_freqs ^ 2)
      N2 <- alleles_out$N2
      Hexp_nei <- (N2 / (N2 - 1)) * het_exp

      genotypes_rand <- genotype_summary2(your_df, k, k + 1)	
      t_f_rand <- genotypes_rand[, 1] == genotypes_rand[, 2]
      homo_obs_rand <- sum(t_f_rand * genotypes_rand[, 3]) / N
      het_obs_rand <- 1 - homo_obs_rand
      f_permuted <- round(1 - (het_obs_rand / Hexp_nei), 3)
      f_vect <- append(f_vect, f_permuted) 
      k = k + 2
    }

    f_avg_vect[j] <- mean(f_vect)
  }

  f_table <- f_multi(your_df)
  f_avg_obs <- unname(f_avg(f_table)[2])
  f_avg_vect <- append(f_avg_vect, f_avg_obs)

  if(f_avg_obs >= 0){
    P_val <- sum(f_avg_vect >= f_avg_obs) / length(f_avg_vect)
  } else if(f_avg_obs < 0){
    P_val <- sum(f_avg_vect < f_avg_obs)
  }

return(c(f_avg = f_avg_obs, P_value = P_val))
}
# End


f_test_avg(geno, 9999)
f_test_avg_out <- f_test_avg(geno, 9999)


write.csv(f_test_avg_out, "f_test_avg.csv") # If you want to output a table to Excel
                                        # Output table will be in your working directory

