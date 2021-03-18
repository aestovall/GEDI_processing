library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

parLapply(cl, file_ls, function(x){
  
  source("R/GEDI_setup.R")
  filename<-gsub(".h5",".txt", 
                 gsub(GEDI.dir,output.dir, x) )
  #if you want to save the GEDI shots into the workspace you can.
  # GEDI.shots<-preprocess.GEDI(input.file=x)
  
  #otherwise save a txt files for bulk processing
  try(preprocess.GEDI(input.file=x, 
                      filename=filename), silent=TRUE)
})

stopCluster(cl)
