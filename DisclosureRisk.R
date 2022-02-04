
# Project : Evaluating disclosure risk and data utility of partially synthetic data by comparing it to local suppression

# Packages needed
library(readr)
library(reghelper)
library(synthpop)
library(sdcMicro)
library(sampling)
library(magrittr)
library(dplyr)
library(caret)
library(nnet)
library(udpipe)
library(Rcpp)
library(tibble)
library(tidyr)
library(fastLink)


# increasing the memory limit
# memory.limit(9999999999)


# Loading the data
adult <- readr::read_csv("C:/Users/mathi/Downloads/adult.data",col_names = FALSE) 

# Converting the data into data frame
obsdata <- as.data.frame(adult)

# Create unique identifier for all observations
obsdata <-obsdata%>%
  udpipe::unique_identifier(fields = colnames(obsdata), start_from = 1L)%>%
  cbind(obsdata)

# renaming the variables
obsdata%<>%
  dplyr::rename(c("unique_id"=".","Age"="X1","workClass"="X2","fnlwgt"="X3",
                  "education"="X4","educationNum"="X5","maritalStatus"="X6",
                  "occupation"="X7","relationship"="X8","race"="X9", "Sex"="X10",
                  "capitalgain"="X11","capitalloss"="X12","hoursPerWeek" = "X13",
                  "nativeCountry"="X14","IncomeClass"="X15"))

#Cleaning the data and removing unrecorded variable values of observations 
obsdata$nativeCountry[which(obsdata$nativeCountry=="?")] <-NA
obsdata$workClass[which(obsdata$workClass=="?")] <-NA
obsdata$occupation[which(obsdata$occupation=="?")] <-NA
obsdata$race[which(obsdata$race=="Other")] <-NA

obsdata%<>%
  na.omit()

# Converting characters into factors
obsdata$workClass <- as.factor(obsdata$workClass)
obsdata$education <- as.factor(obsdata$education)
obsdata$maritalStatus <- as.factor(obsdata$maritalStatus)
obsdata$occupation <-  as.factor(obsdata$occupation)
obsdata$relationship <- as.factor(obsdata$relationship)
obsdata$race <- as.factor(obsdata$race)
obsdata$Sex <- as.factor(obsdata$Sex)
obsdata$nativeCountry <- as.factor(obsdata$nativeCountry)
obsdata$IncomeClass <- as.factor(obsdata$IncomeClass)

dplyr::glimpse(obsdata)
summary(obsdata)

# quasi-identifiers 
quasi_id <- c("Age","Sex","race","maritalStatus","IncomeClass")

# Plotting a histogram of the quasi identifiers
#Age
obsdata%>%
  ggplot2::ggplot(aes(x=Age))+
  ggplot2::geom_histogram(bins = 30, size = 1,  fill="darkgreen", alpha=0.7)

#Sex
obsdata%>%
  ggplot2::ggplot(aes(x= Sex))+
  ggplot2::geom_histogram(stat = "count", size = 1,  fill="darkgreen", alpha=0.7, width=0.4)

#Race
obsdata%>%
  ggplot2::ggplot(aes(x= race))+
  ggplot2::geom_histogram(stat = "count", size = 1,  fill="darkgreen", alpha=0.7, width=0.8)

# Marital status
obsdata%>%
  ggplot2::ggplot(aes(x= maritalStatus))+
  ggplot2::geom_histogram(stat = "count", size = 1,  fill="darkgreen", alpha=0.7, width=0.4)

#Income class
obsdata%>%
  ggplot2::ggplot(aes(x= IncomeClass))+
  ggplot2::geom_histogram(stat = "count", size = 1,  fill="darkgreen", alpha=0.7, width=0.4)


# Monte-Carlo Simulation for the simulation study
# Before we create our function, we check if all the single steps work properly.

n <- 9000
srs_sample <- sampling::srswor(n, nrow(obsdata))
sample_data <- obsdata[srs_sample>0, ]
data <-as.data.frame(sample_data)     # removing the ID column  before applying the methods      

# Frequency table of the variables            
freqCalc(data, keyVars = quasi_id) # checking number of observations violating k-anonymity 

###### Generating 25 different synthetic dataset, thus m=25 and removing the ID column
syn_data <-synthpop::syn(data, m= 25, proper = TRUE,visit.sequence = c(16,2,7,10,11), 
                         drop.not.used = FALSE)

# Comparing synthetic dataset with original dataset
syn_compare<-synthpop::compare(syn_data, data, vars = quasi_id)

# generalized linear model using quasi-identifiers 
syn_mod <- synthpop::glm.synds(IncomeClass ~ Age + Sex + race ,family = "binomial", syn_data)

# Compare synthetic data with original data using the model fits
syn_compare_fit <-synthpop::compare.fit.synds(syn_mod,data , return.plot = F) 

syn_compare_fit


###### Applying local suppression
num_vars = c("Age")

key_vars <- c("Sex","race","maritalStatus","IncomeClass")

obsdata_sdc <- sdcMicro::createSdcObj(data, keyVars = key_vars, numVars = num_vars) 

# sum of unique combinations in data before data is suppressed
sum(get.sdcMicroObj(obsdata_sdc, type = "risk")$individual==1) 

#Suppressing data
res_ls <- obsdata_sdc %>% localSuppression(k=5)

#Percentage of supressions per key variables  
print(res_ls, 'ls')
print(res_ls, "risk")

sum(res_ls@risk$individual[,1])

# data frame of the locally suppressed variables
data_qi <- cbind(res_ls@manipNumVars[,num_vars],res_ls@manipKeyVars[,key_vars])%>%
  dplyr::rename("Age" = "res_ls@manipNumVars[, num_vars]")# renaming the column to Age
colnames(data_qi)

# Replacing the quasi identifiers in the original data with the quasi identifiers suppressed 
data_ls <- res_ls@origData%>%  
  dplyr::select(-c(2,7,10,11,16))%>% #c(1,4,6,9,10,15)
  cbind(data_qi)

# Modeling IncomeClass using the suppressed data
mod_ls <- glm(IncomeClass ~ Age + Sex  + race , family = "binomial", data = data_ls)
summary(mod_ls)
CI.mod_ls <-confint(mod_ls, level = 0.95)
summary(mod_ls)$coef

# Modeling IncomeClass using the original data
mod_origdata <- glm(IncomeClass ~ Age + Sex + race ,family = "binomial", data = res_ls@origData)
summary(mod_origdata)
CI.origdata <- confint(mod_origdata, level = 0.95)

# standardized regression model of local suppressed data
st.mod_ls <-reghelper::beta(mod_ls, x=TRUE, y= FALSE, skip = NULL) 

# standardized regression model of original data
st.mod_origdata <-reghelper::beta(mod_origdata, x=TRUE, y= FALSE, skip = NULL)

# Calculating the confidence interval overlap
CI_overlap <- function(CI_ls, CIobs) {
  # overlap per https://doi.org/10.1198/000313006X124640
  l_over <- max(CI_ls[1], CIobs[1])
  u_over <- min(CI_ls[2], CIobs[2])
  # calculate overlap in units of interval lengths
  Q <- 0.5 * (
    ((u_over - l_over) / (CIobs[2] - CIobs[1])) +
      ((u_over - l_over) / (CI_ls[2] - CI_ls[1]))
  )
  # return results
  return(Q)
}


# create results matrix
Q <- matrix(data = NA,  ncol = nrow(CI.origdata), dimnames=list(NULL, rownames(CI.origdata)))

for (q in 1:nrow(CI.origdata)) {
  Q[1,q]<- CI_overlap(CI_ls = CI.mod_ls[q, ], CIobs = CI.origdata[q, ])
}

# combining the standardized coefficients of the estimates from 
# both local suppressed data and original data and CI overlap
Results_ls <-cbind(coef(mod_ls),coef(mod_origdata))%>%
  cbind(coef(st.mod_ls)[,1] - coef(st.mod_origdata)[,1], t(Q)) %>%
  as.data.frame()%>%
  dplyr::rename(c("local supp"="V1","observed"="V2","Std.Coef_diff"="V3","CI Overlap"="V4"))


# Computing the risk of identification for both synthetic and suppressed data
Obs <- data[,1]%>%cbind(data[quasi_id])%>%rename("unique_id"=".") # Observed data containing only the quasi identifiers

## computing the number of matches using data containing only the quasi identifiers

# Create results matrix for the fraction of matches between observed data and synthesized data
fraction_match_syn <-matrix(data = NA, ncol = 2, dimnames=list(NULL, colnames("Number")))

for (i in 1:2){
  comp_Obs_syn <- syn_data$syn[[i]]%>%
    fastLink::fastLink(Obs, varnames = quasi_id, n.cores = 1)
  
  matched <- getMatches(Obs,syn_data$syn[[i]],comp_Obs_syn,combine.dfs = FALSE) # matches
  dA<-matched$dfA.match; dB<-matched$dfB.match
  commonvars <- intersect(colnames(dA), colnames(dB)) # common variables in both data frames
  colnames(dA)[colnames(dA) %in% commonvars] <-  # Ensuring uniqueness of the column names in observed data
    paste("A.", colnames(dA)[colnames(dA) %in% commonvars], sep="") 
  colnames(dB)[colnames(dB) %in% commonvars] <-  # Ensuring uniqueness of the column names in synthesized data
    paste("B.", colnames(dB)[colnames(dB) %in% commonvars], sep="")
  merged_match <-cbind.data.frame(dA, dB)  # combining both data frames
  mean(merged_match[,"A.unique_id"] != merged_match[,"B.unique_id"]) # Fraction of mismatches
  fraction_match_syn[i] <-1-mean(merged_match[,"A.unique_id"] != merged_match[,"B.unique_id"])# Fraction of matches
}

# Suppressed data containing only the quasi identifiers
Sup_data <- data_ls

# Create results matrix for the fraction of matches between observed data and suppressed data
fraction_match_ls <-matrix(data = NA, dimnames=list(NULL, colnames("Number")))

# comparing the observed data to the suppressed data
comp_Obs_sup <- Sup_data%>%
  fastLink::fastLink(Obs, varnames = quasi_id, n.cores = 1)


matched_ls <- fastLink::getMatches(Obs,Sup_data,comp_Obs_sup,combine.dfs = FALSE) # matches 
dA_obs<-matched_ls$dfA.match; dB_ls<-matched_ls$dfB.match
commonvars_ls <- intersect(colnames(dA_obs), colnames(dB_ls)) # common variables in both data frames
colnames(dA_obs)[colnames(dA_obs) %in% commonvars_ls] <-  # Ensuring uniqueness of the column names in observed data
  paste("A.", colnames(dA_obs)[colnames(dA_obs) %in% commonvars_ls], sep="")
colnames(dB_ls)[colnames(dB_ls) %in% commonvars_ls] <-   # Ensuring uniqueness of the column names in suppressed data
  paste("B.", colnames(dB_ls)[colnames(dB_ls) %in% commonvars_ls], sep="")
merged_ls <-cbind.data.frame(dA_obs, dB_ls) # combining both data frames
mean(merged_ls[,"A.unique_id"] != merged_ls[,"B.unique_id"]) # Fraction of mismatches
fraction_match_ls[i] <-1-mean(merged_ls[,"A.unique_id"] != merged_ls[,"B.unique_id"]) # Fraction of matches





### Now we implement this as a function, which allows to perform our Monte-Carlo-Simulation

# Calculating the confidence interval overlap
CI_overlap <- function(CI_ls, CIobs) {
  # overlap per https://doi.org/10.1198/000313006X124640
  l_over <- max(CI_ls[1], CIobs[1])
  u_over <- min(CI_ls[2], CIobs[2])
  # calculate overlap in units of interval lengths
  Q <- 0.5 * (
    ((u_over - l_over) / (CIobs[2] - CIobs[1])) +
      ((u_over - l_over) / (CI_ls[2] - CI_ls[1]))
  )
  # return results
  return(Q)
}



####  Monte Carlo simulation 

my_simulation <- function( n=n) { 
  
  quasi_id <- c("Age","Sex","race","maritalStatus","IncomeClass")
  
  data <- obsdata[sampling::srswor(n, nrow(obsdata))>0, ] %>% as.data.frame()
  
  syn_data <- synthpop::syn(data, m=25, proper = TRUE, visit.sequence = c(16,2,7,10,11),
                            drop.not.used = FALSE)
  
  # Comparing synthetic dataset with original dataset
  syn_compare<- synthpop::compare(syn_data, data, vars = quasi_id)
  
  # generalized linear model using quasi-identifiers 
  syn_mod <-  synthpop::glm.synds(IncomeClass ~ Age + Sex + race ,family = "binomial", syn_data)
  
  # Compare synthetic data with original data using the model fits
  syn_compare_fit <-syn_mod %>%
    synthpop::compare.fit.synds(data, return.plot = F)
  
  
  ##### Local suppression
  num_vars = c("Age")
  
  key_vars <- c("Sex","race","maritalStatus","IncomeClass")
  
  # creating an sdc object
  data_sdc <- sdcMicro::createSdcObj(data, keyVars = key_vars, numVars = num_vars, weightVar = "fnlwgt")
  
  res_ls <-data_sdc %>% sdcMicro::localSuppression(k=5) # suppression of quasi-identifiers
  
  # data frame of the locally suppressed variables
  data_qi <- cbind(res_ls@manipNumVars[,num_vars],res_ls@manipKeyVars[,key_vars])%>%
    dplyr::rename("Age" = "res_ls@manipNumVars[, num_vars]")# renaming the column to Age
  colnames(data_qi)
  
  # Replacing the quasi identifiers in the original data with the quasi identifiers suppressed 
  data_ls <- res_ls@origData %>%  
    dplyr::select(-c(2,7,10,11,16))%>%
    cbind(data_qi)
  
  # Modeling IncomeClass using the suppressed data
  mod_ls <- glm(IncomeClass ~ Age + Sex + race ,family = "binomial",data = data_ls, control=glm.control(maxit=50))
  CI.mod_ls <-confint(mod_ls, level = 0.95)
  
  # Modeling with the original data
  mod_origdata <- glm(IncomeClass ~ Age + Sex + race ,family = "binomial", data = res_ls@origData, control=glm.control(maxit=50))
  CI.origdata <- confint(mod_origdata, level = 0.95)
  
  # standardized regression model of locally suppressed data
  st.mod_ls <-reghelper::beta(mod_ls, x=TRUE, y= FALSE, skip = NULL) 
  
  # standardized regression model of original data
  st.mod_origdata <-reghelper::beta(mod_origdata, x=TRUE, y= FALSE, skip = NULL)
  
  # create results matrix
  Q <- matrix(data = NA,  ncol = nrow(CI.origdata), dimnames=list(NULL, rownames(CI.origdata)))
  
  
  for (q in 1:nrow(CI.origdata)) {
    Q[1,q]<- CI_overlap(CI_ls = CI.mod_ls[q, ], CIobs = CI.origdata[q, ])
  }
  
  # combining the standardized coefficients of the estimates from both local suppressed data and original data and CI overlap
  Results_ls  <-cbind(coef(mod_ls),coef(mod_origdata))%>%
    cbind(coef(st.mod_ls)[,1] - coef(st.mod_origdata)[,1], t(Q)) %>%
    as.data.frame()%>%
    dplyr::rename(c("local supp"="V1","observed"="V2","Std.Coef_diff"="V3","CI Overlap"="V4"))
  
  
  #### Estimating the risk of identification for both synthetic and suppressed data
  
  Obs <- data # Observed data containing only the quasi identifiers
  
  # Suppressed data containing only the quasi identifiers
  Sup_data <- data_ls[,1]%>%cbind(data_ls[quasi_id])%>%rename("unique_id"=".") 
  
  # compute linked data between observed and suppressed data
  comp_Obs_sup<- fastLink::fastLink(Obs,Sup_data, varnames = quasi_id)
  
  # get matches in both data
  matched_ls <- fastLink::getMatches(Obs,Sup_data,comp_Obs_sup,combine.dfs = FALSE) 
  dA_obs<-matched_ls$dfA.match; dB_ls<-matched_ls$dfB.match
  
  # common variables in both data frames
  commonvars_ls <- intersect(colnames(dA_obs), colnames(dB_ls)) 
  
  # Ensuring uniqueness of the column names in observed data
  colnames(dA_obs)[colnames(dA_obs) %in% commonvars_ls] <-  
    paste("A.", colnames(dA_obs)[colnames(dA_obs) %in% commonvars_ls], sep="")
  
  # Ensuring uniqueness of the column names in suppressed data
  colnames(dB_ls)[colnames(dB_ls) %in% commonvars_ls] <-   
    paste("B.", colnames(dB_ls)[colnames(dB_ls) %in% commonvars_ls], sep="")
  
  # combining both data frames
  merged_ls <-cbind.data.frame(dA_obs, dB_ls) 
  
  # Fraction of mismatches between observed data and suppressed data
  mean(merged_ls[,"A.unique_id"] != merged_ls[,"B.unique_id"]) 
  
  # Fraction of matches observed data and suppressed data
  fraction_match_ls <-1-mean(merged_ls[,"A.unique_id"] != merged_ls[,"B.unique_id"]) 
  
  
  # Create results matrix for the fraction of matches between observed data and synthesized data
  fraction_match_syn <-matrix(data = NA, ncol =syn_data$m, dimnames=list(NULL, colnames("Number")))
  
  # Create results matrix for the mean fraction of matches between observed data and synthesized data
  fraction_mean_match_syn <-matrix(data = NA, ncol = syn_data$m, dimnames=list(NULL, colnames("Number")))
  
  #### compute linked data between observed and synthetic data
  for (i in 1:syn_data$m){
    comp_Obs_syn <- syn_data$syn[[i]]%>%
      fastLink::fastLink(Obs, varnames = quasi_id, n.cores = 1)
    
    # Get matches from both data
    matched <- fastLink::getMatches(Obs,syn_data$syn[[i]],comp_Obs_syn,combine.dfs = FALSE) 
    dA<-matched$dfA.match; dB<-matched$dfB.match
    
    # common variables in both data frames
    commonvars <- intersect(colnames(dA), colnames(dB))
    
    # Ensuring uniqueness of the column names in observed data
    colnames(dA)[colnames(dA) %in% commonvars] <-  
      paste("A.", colnames(dA)[colnames(dA) %in% commonvars], sep="") 
    
    # Ensuring uniqueness of the column names in synthesized data
    colnames(dB)[colnames(dB) %in% commonvars] <-  
      paste("B.", colnames(dB)[colnames(dB) %in% commonvars], sep="")
    
    # combining both data frames
    merged_match <-cbind.data.frame(dA, dB) 
    
    # Fraction of mismatches between observed and synthetic data
    mean(merged_match[,"A.unique_id"] != merged_match[,"B.unique_id"]) 
    
    # Fraction of matches between observed and synthetic data
    fraction_match_syn[i] <-1-mean(merged_match[,"A.unique_id"] != merged_match[,"B.unique_id"])
    fraction_mean_match_syn[i] <- apply(X = as.matrix(fraction_match_syn[i]), MARGIN = 1, FUN = mean)
  }
  
  my_list <- list(syn_compare_fit=syn_compare_fit,Results_ls=Results_ls,fraction_match_ls=fraction_match_ls,
                  fraction_mean_match_syn=fraction_mean_match_syn)
  
  return(my_list)
  
}


set.seed(1001)
result <- list()# Results of each simulation run will be stored here
system.time(for (i in 1:1000){# Number of Simulation runs
  print(i)
  resultX <-my_simulation(n=5000) #With Sample size n
  result[[length(result)+1]] <- resultX
})

# Combining mean results and saving in a data frame
for (i in 1:length(result)){
  #Observed coeff
  Observed <-apply(X = as.matrix(result[[i]]$syn_compare_fit$coef.diff[,2]), MARGIN = 1, FUN = mean)
  
  #Synthetic coeff
  Synthetic <-apply(X = as.matrix(result[[i]]$syn_compare_fit$coef.diff[,1]), MARGIN = 1, FUN = mean)
  
  # Local suppressed coeff
  Local_supp <-apply(X = as.matrix(result[[i]]$Results_ls[,1]), MARGIN = 1, FUN = mean)
  
  #Std.Coef_diff between observed and synthetic
  Std_obs_syn <-apply(X = as.matrix(result[[i]]$syn_compare_fit$coef.diff[,4]), MARGIN = 1, FUN = mean)
  
  # Std.Coef_diff between observed and Local suppressed
  Std_obs_Ls <-apply(X = as.matrix(result[[i]]$Results_ls[,3]), MARGIN = 1, FUN = mean)
  
  #Synthetic - confidence interval overlap
  ci_syn <-apply(X = as.matrix(result[[i]]$syn_compare_fit$ci.overlap), MARGIN = 1, FUN = mean) 
  
  # Local suppressd - confidence interval overlap
  ci_Ls <-apply(X = as.matrix(result[[i]]$Results_ls[,4]), MARGIN = 1, FUN = mean)
  
  #fraction of matches in suppressed data
  # mean of the fraction of matches in suppressed data
  mean_fraction_match_ls <- apply(X = as.matrix(result[[i]]$fraction_match_ls),MARGIN = 2, FUN = mean)
  
  #fraction of matches in synthetic data
  mean_match_syn <- apply(X = as.matrix(result[[i]]$fraction_mean_match_syn),MARGIN = 2, FUN = mean) 
  # mean of the fraction of matches in synthetic data
  mean_fraction_match_syn <- apply(X = as.matrix(mean_match_syn),MARGIN = 2, FUN = mean)
  
  Fin_res <- data.frame(Observed,Synthetic,Local_supp,Std_obs_syn,Std_obs_Ls,ci_syn,ci_Ls,mean_fraction_match_ls,mean_fraction_match_syn)
  
}


Fin_res # data frame of combined results

mean(Fin_res[,6]) # Total mean confidence interval overlap of synthetic data
mean(Fin_res[,7]) # Total mean confidence interval overlap of suppressed data


# Save results in csv file format
write.csv(Fin_res,"C:\\Users\\mathi\\Documents\\sim_res.csv") #combined results in a data frame



# A line plot of the confidence interval overlap estimates
variable_names <- c("Intercept","Age","SexMale","raceAsian-Pac-Islander","raceBlack","raceWhite") # Key_variables
new_data <-Fin_res%>%
  tibble::add_column(variable_names)%>%
  tidyr::pivot_longer(cols = c(ci_syn,ci_Ls),values_to = "Confidence interval overlap", names_to = "Method")
new_data%>%
  ggplot2::ggplot(aes(y=`Confidence interval overlap`, x=variable_names,color=Method)) +
  geom_point(alpha = 1,size=1.5)+
  geom_line(aes(colour = Method, group = Method), size=0.5)+
  expand_limits(y=0)



# A line plot of the standard regression coefficient estimates
variable_names <- c("Intercept","Age","SexMale","raceAsian-Pac-Islander","raceBlack","raceWhite") # Key_variables
new_data <-Fin_res%>%
  tibble::add_column(variable_names)%>%
  tidyr::pivot_longer(cols = c(Observed,Synthetic,Local_supp),values_to = "standard regression coefficients", names_to = "Method")
new_data%>%
  ggplot2::ggplot(aes(y=`standard regression coefficients`, x=variable_names,color=Method)) +
  geom_point(alpha = 1)+
  geom_line(aes(colour = Method, group = Method), size=0.5)+
  expand_limits(y=0)







