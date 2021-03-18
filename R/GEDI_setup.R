library(h5)
source("R/GEDI_processing_FUN.R")
beams<-c("0000","0001","0010","0011",
         "0101","0110", "1000", "1011")

GEDI.dir<-"D:/GEDI/"
output.dir<-"D:/GEDI/output/"

file_ls<-subset.files(GEDI.dir, output.dir)
