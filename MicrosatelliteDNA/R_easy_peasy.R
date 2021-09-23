### Run this code from start to end.
### It will load all the required functions and output the two comma delimited...
### ...files into your working directory that you need for the homework.

### The output comma delimited files can be opened in MS Excel.
### Don't forget to interpret the results (either typed or written).
### Provide the annotated printout when you turn in your homework.


### REMEMBER to put "Sparrow_only.csv" into your working direcotry.
### Your working directory is typically your Documents folder.

###### Start #######

geno <- read.csv("Sparrow_only.csv")

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


f_avg <- function(fs){
	mean_f_raw <- mean(fs$f_raw)
	mean_f_corr <- mean(fs$f_corrected)
	
return(c("mean f_raw"= mean_f_raw, "mean f_corr" = mean_f_corr))
}


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

f_tests <- f_test_multi(geno, 9999)
write.csv(f_tests, "f_tests.csv")

f_test_avg_out <- f_test_avg(geno, 9999)
write.csv(f_test_avg_out, "f_test_avg.csv") # If you want to output a table to Excel
   
###### End #######







