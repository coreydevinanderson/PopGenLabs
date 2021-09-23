### Running the code in this file will import the uDNA data into R. ###

## uDNA data must be in diploid two column format where the first column...
## ...contains IDs for each individual that was sampled.

## There is code to import either the armadillo or sparrow data.
## There is also a function at the end that will provide information about...
## ...each file.

## If you are using a .txt file, you must copy the code and paste it into...
## ...the R Editor. To open the R Editor box, select the "File" tab, then...
## ..."New script", then paste your code into the box.

## If you are using an .R file, you must open it in R by selecting "File"...
## ...then "Open script".

## You can also use ctrl+N (for "New script") and ctrl+O ("Open script").

## To run a block of code, highlight it, then hit: ctrl+R
## You can also highlight the code and then choose: "Run line or selection"...
## ..from the "Edit" tab.

## For Macintosh, use command+Enter to run a block of selected code.

##########

# To import the data, the following files must be in your working directory:

# 1) Armadillo_only.csv
# 2) Sparrow_only.csv

# To determine your working directory:

getwd()

##########

rm(list = ls()) # Clear R's memory


##################
# Armadillo data #
##################

# Run this code (from # Start to # End) to import "Armadillo_only.csv"

# Start
geno_raw <- read.csv("Armadillo_only.csv")  #read in the file
geno_adult <- subset(geno_raw, Age.class == "Adult")  #select only adults
geno <- geno_adult[, - dim(geno_adult)[2]]  #remove the last column
# End


########################
# Grasshopper sparrows #
########################

# Run this code (from # Start to # End) to import "Sparrow_only.csv"

# Start
geno <- read.csv("Sparrow_only.csv")
# End


################
## summary_df ##
################

## Description: provides a simple summary of your imported genotype file.

## Arguments: 
# your_df: a data.frame with your genotype data (= geno)

# Run the code (from # Start to # End and then call the function)

# Start
summary_df <- function(your_df){

  df_dim <- dim(your_df)
  
  line1 <- paste("Your data.frame has ", df_dim[1], " rows.", sep = "")
  line2 <- "The data are in diploid two-column format: two columns per locus."
  line3 <- paste("There are ", (df_dim[2] - 1 )/ 2, " loci; the first locus is in columns 2 & 3.", sep = "")
  line4 <- paste("The last locus is in columns is ", df_dim[2] - 1, " & ", df_dim[2], ".", sep = "")
 
  lines <- c(line1, line2, line3, line4)

return(cat(lines, sep = "\n"))
}
# End

summary_df(geno) # Call the function and give it "geno" as an argument


