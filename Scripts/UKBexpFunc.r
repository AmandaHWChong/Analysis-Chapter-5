#R function to extract exposures of interest from UK Biobank ICD10 data derived from script created by Tom Dudding and Kimberley Burrows

#As the observational analyses uses self-reported, ICD9, and ICD10 data to identify cases for exposure of interest, the analyses could be subject to reverse causation. Therefore sensitivity analyses will be performed by identifying cases after initial platelet phenotype measurements (platelet count and mean platelet volume) to identify any evidence of reverse causation. 
#Exposures of interest: Acute myocardial infarction, deep vein thrombosis, pulmonary embolism 

# Options
# dat - UK Biobank dataframe - the variable names should ideally match those provided by UK Biobank, NB the id column must be projectID
# exposureCode - vector of ICD codes of interest, for ICD10 this is a vector of class character.
# exp_col_start - the first column including the exposure data (the default corresponds to the default name of the ICD10 data)
# exp_col_fin - the final column including the exposure data.
# date_col_start - first column of the ICD10 date variable
# date_col_fin - last column of the ICD10 date variable

UKBexpFunc<-function(dat,exposureCode, sitename, other=F, 
                        exp_col_start= f.41270.0.0,
                        exp_col_fin=f.41270.0.212,
                        date_col_start = f.41280.0.0,
                        date_col_fin =f.41280.0.225) {
  
  if (other ==T) {
    #generate list of ICDcodes not inlcuding those of interest
    allICD<-unlist(dplyr::select(dat,get(exp_col_start):get(exp_col_fin)), use.names = F)
    allICD<-allICD[!is.na(allICD)]
    allICD<-allICD[!duplicated(allICD)]
    exposureCode<-allICD[!allICD %in% exposureCode]
  }
  
  #idenitfy all participants with specific ICD10 codes
  exposure_cases<-dplyr::select(dat,projectID,get(exp_col_start):get(exp_col_fin))
  exposureDate<-dplyr::select(dat,projectID,get(date_col_start):get(date_col_fin))
  #rownames(exposure_cases)<-rownames(dat)
  #rownames(exposureDate)<-rownames(dat)
  exposures2<-exposure_cases[apply(exposure_cases, 1, function(r) any(r %in% exposureCode)),]
  if (length(exposures2$projectID)<1) {
    message("No exposure of this type")
    return(dat)
    } else {
    
  exposureDate2<-exposureDate[apply(exposure_cases, 1, function(r) any(r %in% exposureCode)),]
    
  #this loop takes the exposure codes and dates and for each participant generates a list of dates that they were diagnosed with the exposure of interest
  exposureList<-as.list(NULL)
  for (i in 1:length(exposures2$projectID)) {
    person1C<-unlist(exposures2[i,-1])
    person1D<-as.Date(unlist(exposureDate2[i,-1]),origin = "1992-03-31") #also converts to r internal dates (number of days since 1992-03-31)
    #set to missing non interest exposure codes
    person1C<-ifelse(person1C %in% exposureCode,person1C,NA) #keeps the NAs
    #person1C<-person1C[person1C %in% oral] #removes the NAs
    #remove dates for exposure codes not interested in
    person1D<-as.Date(ifelse(!is.na(person1C),person1D,NA),origin = "1992-03-31")
    person1D<-as.Date(person1D[!is.na(person1C)],origin = "1992-03-31")
     
    #remove duplicated exposures
    if (length(person1D)>1) {
      person1C<-person1C[!duplicated(person1D)]
      person1D<-as.Date(person1D[!duplicated(person1D)],origin = "1992-03-31")
      #sort by earliest diagnosis first
      iD <- order(person1D)
      person1D<-as.Date(person1D[iD],origin = "1992-03-31")
      person1C<-person1C[iD]
            
    } else {
      person1D<-as.Date(person1D,origin = "1992-03-31")
    }
    exposureList[[i]]<-paste(person1C,person1D, sep  = "/")
    #exposureList[[i]]<-c(person1C,as.Date(person1D,origin = "1992-03-31")
  }
  
  #create a dataframe of exposure diagnosis dates
  
  exposureDF<-plyr::ldply(exposureList, rbind)
  maxC<-max(unlist(lapply(exposureList,length)))
  #exposureDF[,1:maxC]<-as.character(exposureDF[,1:maxC])
  colnames(exposureDF)<-paste0(sitename,1:maxC)

  
  exposureDF$projectID<-exposure2$projectID
  
  
  
  
  #merge exposure data back into fulldataset
  fulldat<-merge(dat,exposureDF, by="projectID", all.x = T)
  return(fulldat)
    }
}
