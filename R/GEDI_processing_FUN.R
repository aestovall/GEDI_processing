preprocess.GEDI<-function(input.file,
                          beams=c("0000","0001","0010","0011",
                                  "0101","0110", "1000", "1011"),
                          filename=NA,
                          ....){
  
  file <- try(h5::h5file(name = input.file,
                         mode = "a"), silent = TRUE)
  
  if(class(file)=="H5File"){
    data.sets<-list.datasets(file)
    beams_ls<-lapply(1:8, function(i){
      require(h5)
      lat<-file[paste("/BEAM",beams[i],"/geolocation/lat_highestreturn", sep="")]
      if(is.null(nrow(lat))) next
      lon<-file[paste("/BEAM",beams[i],"/geolocation/lon_highestreturn", sep="")]
      elev_lowestmode<-file[paste("/BEAM",beams[i],"/geolocation/elev_lowestmode", sep="")]
      elev_TDX<-file[paste("/BEAM",beams[i],"/geolocation/digital_elevation_model", sep="")]
      
      
      rh100<-file[paste("/BEAM",beams[i],"/rh100", sep="")]
      cover<-file[paste("/BEAM",beams[i],"/cover", sep="")]
      pai<-file[paste("/BEAM",beams[i],"/pai", sep="")]
      fhd<-file[paste("/BEAM",beams[i],"/fhd_normal", sep="")]
      num_detectedmodes<-file[paste("/BEAM",beams[i],"/num_detectedmodes", sep="")]
      modis_treecover<-file[paste("/BEAM",beams[i],"/land_cover_data/modis_treecover", sep="")]
      
      
      q1<-file[paste("/BEAM",beams[i],"/l2a_quality_flag", sep="")]
      surface_flag<-file[paste("/BEAM",beams[i],"/surface_flag", sep="")]
      sensitivity<-file[paste("/BEAM",beams[i],"/sensitivity", sep="")]
      shot_number<-file[paste("/BEAM",beams[i],"/shot_number", sep="")]
      
      delta_time<-file[paste("/BEAM",beams[i],"/geolocation/delta_time", sep="")]
      
      return(data.frame(lon=lon[1:nrow(lon)],
                        lat=lat[1:nrow(lat)],
                        elev_lowestmode=elev_lowestmode[1:nrow(elev_lowestmode)],
                        elev_TDX = elev_TDX[1:nrow(elev_TDX)],
                        rh100=rh100[1:nrow(rh100)]/100,
                        cover=cover[1:nrow(cover)],
                        pai=pai[1:nrow(pai)],
                        fhd=fhd[1:nrow(fhd)],
                        num_detectedmodes=num_detectedmodes[1:nrow(num_detectedmodes)],
                        modis_treecover=modis_treecover[1:nrow(modis_treecover)],
                        q1=q1[1:nrow(q1)],
                        surface_flag=surface_flag[1:nrow(surface_flag)],
                        sensitivity=round(sensitivity[1:nrow(sensitivity)],1),
                        shot_number=shot_number[1:nrow(shot_number)],
                        delta_time=delta_time[1:nrow(delta_time)],
                        beam = as.character(beams[i]))
             # file = file_ls[ii]
      )
      
    })
    
    if(!is.na(filename)) data.table::fwrite(do.call(rbind,beams_ls), 
                       file = filename, 
                       sep = " ", na = "NA")
    
    return(do.call(rbind,beams_ls))
    
    print(Sys.time()-time_check)
    
  }   else warning("Not h5 file!")
  
}

subset.files<-function(GEDI.dir, output.dir){
  GEDI.files<-list.files(GEDI.dir, pattern = ".h5", full.names = TRUE)
  GEDI.files.name<-list.files(GEDI.dir, pattern = ".h5", full.names = FALSE)
  GEDI.files.output.name<-list.files(output.dir, pattern = ".txt", full.names = FALSE)
  file_ls<-GEDI.files[!(GEDI.files.name %in%
                       gsub(".txt",".h5",GEDI.files.output.name))]
  return(file_ls)
}
