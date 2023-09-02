library(readr)
library(readxl)
library(openxlsx)
library(dplyr)
### necessary libraries for the program

dfmethsite = read_tsv("Suplementary_Data.TSV")
### build a data frame for the methylation sites used

names_dfmethsite = dfmethsite[-c(2:nrow(dfmethsite)), ]
### get rid of all other rows, and only examine the patient/samples in one list - use this for testing purposes

healthycolondf = list()
healthystomachdf = list()
healthylungdf = list()
healthybreastdf = list()
cancercolondf = list()
cancerstomachdf = list()
cancerlungdf = list()
cancerbreastdf = list()
### create lists- use this for testing purposes

for (i in colnames(names_dfmethsite))
  {
  if (startsWith(i, "BC_CRC_ND"))
  {
    healthycolondf[[length(healthycolondf)+1]] = i
  }
  else if(startsWith(i, "BC_Lung_ND"))
  {
    healthylungdf[[length(healthylungdf)+1]] = i
  }
  else if(startsWith(i, "BC_Stomach_ND"))
  {
    healthystomachdf[[length(healthystomachdf)+1]] = i
  }
  else if(startsWith(i, "BC_Breast_ND"))
  {
    healthybreastdf[[length(healthybreastdf)+1]] = i
  }
  else if(startsWith(i, "BC_CRC_TD"))
  {
    cancercolondf[[length(cancercolondf)+1]] = i
  }
  else if(startsWith(i, "BC_Lung_TD"))
  {
    cancerlungdf[[length(cancerlungdf)+1]] = i
  }
  else if(startsWith(i, "BC_Stomach_TD"))
  {
    cancerstomachdf[[length(cancerstomachdf)+1]] = i
  }
  else if(startsWith(i, "BC_Breast_TD"))
  {
    cancerbreastdf[[length(cancerbreastdf)+1]] = i
  }
}
### run through column names and sort each cell into lists- use this for testing purposes

print(healthycolondf)
print(healthylungdf)
print(healthystomachdf)
print(healthybreastdf)
print(cancercolondf)
print(cancerlungdf)
print(cancerstomachdf)
print(cancerbreastdf)
#### print out all the lists generated- use this for testing purposes

healthy_sample_list = list(healthycolondf, healthylungdf, healthystomachdf, healthybreastdf)
cancerous_sample_list = list(cancercolondf, cancerlungdf, cancerstomachdf, cancerbreastdf)
### create big healthy and big cancerous list -> for testing reference

nrows = nrow(dfmethsite)
ncols = ncol(dfmethsite)
cellremoved = 0
### count the cells removed
rownumber = 1
### count the rows 
stop = FALSE
### Set all values; stop value- prevents going into blank cells in excel

for (i in 1:nrows)
{
  if (stop == TRUE)
  {
    break
  }
  rownumber = rownumber + 1
  print(paste("Row Number:", rownumber))
    
    for (j in 1:ncols)
    {
      if (rownumber > 608)
      {
        stop = TRUE
        break()
      }
      
      if (is.na(dfmethsite[i,j]))
      {
        cellremoved = cellremoved + 1
        rownumber = rownumber + 1
        dfmethsite = dfmethsite %>% slice(-c(i))
        print(paste("rows removed equals", cellremoved))
      }
    }
}
### Go through each row, then by row go through columns and check whether cells are empty. If the cells are empty, or null, then remove the entire row. Keep count of rows removed.

write.xlsx(dfmethsite, file = "T-Test_RealMethylation_Values.xlsx")
print("T-test_values have been written")
### Create an excel file with the updated rows and save it

healthydf = select(dfmethsite, -matches("TD"))
healthydf = select(healthydf, -matches("TSH"))
cancerdf = select(dfmethsite, -matches("ND"))
cancerdf = select(cancerdf, -matches("TSH"))
### create a cancerous and healthy data frame based on the naming scheme

t_test_results = vector("list", nrow(healthydf))
### create a vector to hold all the p values of the t-tests

for (i in 1:nrow(healthydf))
{
  t_test_results[[i]] = t.test(healthydf[i,], cancerdf[i,])$p.value
}
### append to the vector the p values for each row, compared between cancerous and healthy data sets

dfmethsite$p_value = t_test_results
### append a column to the main data spreadsheet with the chromosomes that has the p values

dfmethsite = arrange(dfmethsite, desc(p_value))
### sort the p values from greatest to smallest

write.xlsx(dfmethsite, "Complete Data Sheet with P-values.xlsx", overwrite = TRUE)
### write this file

q = 0.30
### FDR

n = length(t_test_resultsbreast)
### Ranking of p-values (updated in the loop)

m = length(t_test_resultsbreast)
### Total number of values 

criticalvalue=0

for (i in dfmethsite$p_value)
{
  nqm = q * (n/m)
  print(nqm)
  print(i)
  
  if (i < nqm)
  {
    criticalvalue = i
    print("done")
    break()
  }
  
  n = n - 1
}
### according to the correction, loop through and find the greatest p-value that is also less than its nqm. Set that as the critical value

dfmethsite = dfmethsite[-c(which(dfmethsite$p_value > criticalvalue)), ]
### delete all rows that have p-values higher than the critical value-> removing all insignificant chromosomal identification

write.xlsx(dfmethsite, "updated_values.xlsx")
### update the sheet with the new values- use this for the regression model

  
