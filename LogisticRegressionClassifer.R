library(readr)
library(readxl)
library(openxlsx)
library(dplyr)
library(tibble)
library(glmnet)
### necessary libraries for the program

file = read_excel("Training_Set_1.xlsx")
### open the file

file = select(file, -c(2:5, 8:ncol(file)))
### delete the unnecessary columns on the total file set- keep only patient ID, Sample Type, and Cancer Tissue

healthyfile = file[file$`Sample Type` == "Healthy", ]
preesophagusfile = file[file$`Cancer Tissue` == "Esophagus" & file$`Sample Type` == "Pre-Diagnosis", ]
prestomachfile = file[file$`Cancer Tissue` == "Stomach" & file$`Sample Type` == "Pre-Diagnosis", ]
precolonfile = file[file$`Cancer Tissue` == "Colon" & file$`Sample Type` == "Pre-Diagnosis", ]
preliverfile = file[file$`Cancer Tissue` == "Liver" & file$`Sample Type` == "Pre-Diagnosis", ]
prelungfile = file[file$`Cancer Tissue` == "Lung" & file$`Sample Type` == "Pre-Diagnosis", ]
postesophagusfile = file[file$`Cancer Tissue` == "Esophagus" & file$`Sample Type` == "Post-Diagnosis", ]
poststomachfile = file[file$`Cancer Tissue` == "Stomach" & file$`Sample Type` == "Post-Diagnosis", ]
postcolonfile = file[file$`Cancer Tissue` == "Colon" & file$`Sample Type` == "Post-Diagnosis", ]
postliverfile = file[file$`Cancer Tissue` == "Liver" & file$`Sample Type` == "Post-Diagnosis", ]
postlungfile = file[file$`Cancer Tissue` == "Lung" & file$`Sample Type` == "Post-Diagnosis", ]
### create data frames for healthy, and types of cancer (pre and post)

dataframelist = list(healthyfile, preesophagusfile, prestomachfile, precolonfile, preliverfile, prelungfile, postesophagusfile, postliverfile, poststomachfile, postlungfile, postcolonfile)
### create a list of all the data frames to run through

df1 = list()
df2 = list()
### create lists to store the chosen and not chosen patients

totaccuracy = 0
ntotaccuracy = 0
counter = 0
### set the total accuracy of chosen, not chosen and counter to zero

while (counter < 100)
{
  
  for(i in dataframelist) 
  {
    i$newcol <- sample(2, nrow(i), replace = TRUE)
    ### randomly select number, and add a new column
    
    i = mutate(i, Chosen=ifelse(newcol==1, "Chosen", "Notchosen"))
    ### based on the randomly generated, either chose for training, or don't
    
    i = select(i, -c(newcol))
    ### delete the newcol row since we're done with it
    
    j = subset(i, Chosen == "Notchosen")
    ### create a sub dataframe of the not chosen
    
    i = subset(i, Chosen == "Chosen") 
    ### create a sub dataframe of the chosen
    
    df1[[length(df1) + 1]] = j
    ### add the not chosen df to the list
    
    df2[[length(df2) + 1]] = i
    ### add the chosen data frame to the list
  }
  ### Go through each data frame, assign its cells random numbers in a new column. Then if it's 1- it's selected for training, if it's 2- it's not selected
  
  df3 = Reduce(function(x, y) merge(x, y, all=TRUE), df2)
  ### Chosen - Combine all the data frames into one big data frame to use for the testing and training set
  
  methylvaluesleavein = Reduce(function(x, y) merge(x, y, all=TRUE), df1)
  ### The not chosen or leave in set - Combine all the data frames into one big data frame
  
  chosesets = sample(c(TRUE, FALSE), nrow(df3), replace=TRUE, prob=c(0.5,0.5))
  trainingset = df3[chosesets,]
  testingset = df3[!chosesets,]
  ### randomly generate a training set and a testing set- equal odds
  
  healthyfile$`Patient ID` = gsub("-", ".", as.character(healthyfile$`Patient ID`))
  trainingset$`Patient ID` = gsub("-", ".", as.character(trainingset$`Patient ID`))
  testingset$`Patient ID` = gsub("-", ".", as.character(testingset$`Patient ID`))
  methylvaluesleavein$`Patient ID` = gsub("-", ".", as.character(methylvaluesleavein$`Patient ID`))
  ### due to poor file formatting- change the hyphens to periods in all the patient ID's for later reference
  
  methylvaluestrain = read_excel("updated_values.xlsx")
  ### taken from the t-tests conducted on the methylation values- has only the significant chromosomes-> use this for training
  
  methylvaluestest = read_excel("updated_values.xlsx")
  ### taken from the t-tests conducted on the methylation values- has only the significant chromosomes-> use this for testing
  
  methylvaluesleave = read_excel("updated_values.xlsx")
  ### taken from the t-tests conducted on the methylation values- has only the significant chromosomes-> use this for the leave in set
  
  for (i in colnames(methylvaluestrain))
  {
    
    if (i %in% testingset$`Patient ID`)
    {
      methylvaluestrain = select(methylvaluestrain, -contains(i))
    }
    ### remove the patient ids that are in the testing set
    
    if (!i %in% trainingset$`Patient ID`& !i%in% testingset$`Patient ID`)
    {
      methylvaluestrain = select(methylvaluestrain, -contains(i))
    }
    ### remove the patient ids that are not in the training, or testing set (ie. the ones not chosen)
    
  }
  ### for training- get the data frame with only the training patient IDs and data
  
  for (i in colnames(methylvaluestest))
  {
    
    if (i %in% trainingset$`Patient ID`)
    {
      methylvaluestest = select(methylvaluestest, -contains(i))
    }
    ### remove the patient ids that are in the testing set
    
    if (!i %in% trainingset$`Patient ID`& !i%in% testingset$`Patient ID`)
    {
      methylvaluestest = select(methylvaluestest, -contains(i))
    }
    ### remove the patient ids that are not in the training, or testing set (ie. the ones not chosen)
    
  }
  ### for testing - get the data frame with only the testing patient IDs and data
  
  for (i in colnames(methylvaluesleave))
  {
    
    if (!i %in% methylvaluesleavein$`Patient ID`)
    {
      methylvaluesleave = select(methylvaluesleave, -contains(i))
    }
    ### remove the patient ids that are not in the leave in set
    
  }
  ### for the leave in set- get the data frame with only the leave in patient IDs and data
  
  methylvalues = as.data.frame(t(methylvaluestrain))
  ### flip the rows and columns by transposition- for training
  
  methylvalues2 = as.data.frame(t(methylvaluestest))
  ### flip the rows and columns by transposition- for testing
  
  methylvalues3 = as.data.frame(t(methylvaluesleave))
  ### flip the rows and columns by transposition- for the leave in set
  
  methylvalues$patientID = colnames(methylvaluestrain)
  ### add a column in for patient id's - use this to tell the algorithm whether the sample is cancerous is not.
  
  methylvalues$Cancer = with(methylvalues, ifelse(!methylvalues$patientID %in% healthyfile$`Patient ID`, "Yes", "No"))
  ### add a column that indicates whether the sample is positive for cancer or not, based on if its patient ID is not in the healthyfile df
  
  methylvalues$nCancer[methylvalues$Cancer == "Yes"] <- 1
  methylvalues$nCancer[methylvalues$Cancer == "No"] <- 0
  ### assign binary values to yes and no for the regression
  
  x <- as.matrix(methylvalues[,2:477])
  z <- as.matrix(methylvalues2[,2:477])
  h <- as.matrix(methylvalues3[,2:477])
  ### matrices of the chromosomal methylation values
  
  y <- as.matrix(methylvalues$nCancer)
  ### matrix of the binary values for cancer (0 or 1)
  
  cv.lasso = cv.glmnet(x, y, alpha = .1, family = "binomial")
  ### calculate the optimal lambda (fitting value for the regression) to use in the regression- mimizing the cross validation error
  
  lsm <- glmnet(x, y, family = binomial, alpha = .1 , lambda = cv.lasso$lambda.min)
  ### create the regression model that compares cancer to all chromosomes sets- uses a Lasso regression, putting less weight on the insignifcant chromosomal markers
  
  testprediction = predict(lsm, z, type="link")
  ### make a prediction for the test set based on the lsm regression model 
  
  leaveoutprediction = predict(lsm, h, type="link")
  ### make a prediction for the leave out set based on the lsm regression model
  
  methylvalues2$prediction = testprediction
  methylvalues3$prediction = leaveoutprediction
  ### append the prediction probabilities to the test set and leave out set data frames
  
  methylvalues2$patientID = colnames(methylvaluestest)
  methylvalues3$patientID = colnames(methylvaluesleave)
  ### add a column in for patient id's - for the test and also leave out data frames- check our accuracy
  
  correct = 0
  incorrect = 0
  ### set correct values for the test set
  
  ncorrect = 0
  nincorrect = 0
  ### set correct values for leave out set
  
  cancerprediction = subset(methylvaluese, methylvaluese$prediction > .565)
  healthyprediction = subset(methylvaluese, methylvaluese$prediction< .565)
  ### create two data frames for the test set- one with what it has labeled to be cancer, one as healthy- cut off point is a probability of .565 (optimal value)
  
  ncancerprediction = subset(methylvalues3, methylvalues3$prediction > .565)
  nhealthyprediction = subset(methylvalues3, methylvalues3$prediction < .565)
  ### create two data frames for the leave out set- one with what it has labeled to be cancer, one as healthy- cut off point is a probability of .565 (optimal value)
  
  for (i in cancerprediction$patientID)
  {
    if (i %in% healthyfile$`Patient ID`)
    {
      incorrect = incorrect + 1
    }
    
    else
    {
      correct = correct + 1
    }
  }
  ### count the number of correct and incorrect cancer diagnoses in the test set results
  
  for (i in healthyprediction$patientID)
  {
    if (i %in% healthyfile$`Patient ID`)
    {
      correct = correct + 1
    }
    
    else
    {
      incorrect = incorrect + 1
    }
  }
  ### count the number of correct and incorrect non-cancer diagnoses in the test set results
  
  for (i in ncancerprediction$patientID)
  {
    if (i %in% healthyfile$`Patient ID`)
    {
      nincorrect = nincorrect + 1
    }
    
    else
    {
      ncorrect = ncorrect + 1
    }
  }
  ### count the number of correct and incorrect cancer diagnoses in the leave out set results
  
  for (i in nhealthyprediction$patientID)
  {
    if (i %in% healthyfile$`Patient ID`)
    {
      ncorrect = ncorrect + 1
    }
    
    else
    {
      nincorrect = nincorrect + 1
    }
  }
  ### count the number of correct and incorrect non-cancer diagnoses in the leave out set results
  
  accuracy = correct / (incorrect + correct)
  naccuracy = ncorrect / (nincorrect + ncorrect)
  ### compute accuracy for the test set and leave out set
  
  totaccuracy = totaccuracy + accuracy
  ntotaccuracy = ntotaccuracy + naccuracy
  ### if looping through total iterations, keep track of total accuracy sum, which will be later referenced.
  
  write.xlsx(methylvalues, "Useforlogisitcregression.xlsx")
  ### write the training set data frame and save it
  
  write.xlsx(methylvaluese, "TestedSet.xlsx")
  ### write the tested set data frame and save it
  
  write.xlsx(methylvalues3, "TestedLeaveinSet.xlsx")
  ### write the leave in set data frame and save it
  
  counter = counter + 1
  ### if you are looping, add to the counter
  
  print(accuracy)
  print(naccuracy)
  print(counter)
  ### if you are looping, keep track of the accuracy per trial, for the test and leave out set, and the counter
  
}

finalaccuracy = totaccuracy/counter
nfinalaccuracy = ntotaccuracy/counter
### compute the final accuracy for the tested set and leave out set

print(paste("The final accuracy for the test set is", finalaccuracy))
print(paste("The final accuracy for the leave out set is", nfinalaccuracy))
### print the final accuracy values

print("Analysis has completed")


