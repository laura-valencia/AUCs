library(tidyverse)
library(ggplot2)
library(stringr)
library(DescTools)  #for AUC

# Data to use: Survival fractions from clonogenic assays. Each treatment has 3 points, forming a curve.


Get_AUCs <- function(rep1, rep2, rep3){
  
  AUC_function <- function(Y_vector){
    
    Result <- AUC(x = c(0,2,3), y = Y_vector)
    
    return(Result)
    
  }
  
  rep1$AUC <- apply(rep1[,-(1:3)], 1, AUC_function)
  rep2$AUC <- apply(rep2[,-(1:3)], 1, AUC_function)
  rep3$AUC <- apply(rep3[,-(1:3)], 1, AUC_function)
  
  Results <- data.frame("Drug" = rep1$Drug,
                        "CellType" = rep1$CellType,
                        "Concentration" = rep1$Concentration,
                        "AUC1" = rep1$AUC,
                        "AUC2" = rep2$AUC,
                        "AUC3" = rep3$AUC)
  
  return(Results)
  
}

main_fun <- function(data1, data2, data3){
  
  
  AUC_DMSO <- Get_AUCs(data1, data2, data3) %>% mutate(Normalization = "DMSO")
  
  AUC_Treat <- Get_AUCs(data1, data2, data3) %>% mutate(Normalization = "Treatment")
  
  Results <- bind_rows(AUC_DMSO, AUC_Treat)
  
  
}

## This part was later prepared into functions to process multiple files and cell lines


# Multiple comparisons for AUCs

# prepare table
Norm_type <- "DMSO"
data <- AUCs_ARHi %>% filter(Normalization==Norm_type) #Select the values normalized by DMSO

# concatenate drug and concentration columns to obtain groups

data$group <- paste(data$Drug, data$Concentration, sep="_")

# gather table

data_g <- gather(data, "AUC1", "AUC2", "AUC3", key="repeat", value = "AUC")

res.aov <- aov(AUC ~ group, data = data_g)
summary(res.aov)

Results <- TukeyHSD(res.aov)
capture.output(Results, file = "anova results.txt")

## comparing between cell lines
Norm_type <- "DMSO"
nSamples <- 24

cell1 <- AUCs_ARHi %>% filter(Normalization==Norm_type)
cell2 <- AUCs_PCDNA %>% filter(Normalization==Norm_type)


cell1$group <- paste(cell1$Drug, cell1$Concentration, cell1$CellType, sep="_")
cell2$group <- paste(cell2$Drug, cell2$Concentration, cell2$CellType, sep="_")

data1_g <- gather(cell1, "AUC1", "AUC2", "AUC3", key="repeat", value = "AUC")
data2_g <- gather(cell2, "AUC1", "AUC2", "AUC3", key="repeat", value = "AUC")

whole_list <- c()


for (i in seq(from=1, to=nSamples, by=3)){
  
  a=i
  b=a+2
  
  
  test_data <- rbind(data1_g[a:b,], data2_g[a:b,])
  
  res.aov <- aov(AUC ~ group, data = test_data)
  results <- summary(res.aov)
  
  cell1_name <- data1_g[a,5]
  cell2_name <- data2_g[a,5]
  print(paste(cell1_name, " VS ", cell2_name))
  print(results)
  
  #save all results that are usually printed to screen into a list of lists.
  Results_element <- c(paste(cell1_name, " VS ", cell2_name), results)
  whole_list <- c(whole_list, Results_element)

}

# Print the whole list of results to a txt file.
sink(paste(data1_g[1,2], "VS", data2_g[1,2], ".txt", sep="_"))
print(whole_list)
sink()
