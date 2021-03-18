

GEDI.processed.files<-list.files(output.dir, pattern="txt", full.names = TRUE)

GEDI<-data.table::fread(GEDI.processed.files[1])

GEDI<-GEDI[GEDI$q1==1,]

library(ggplot2)
library(reshape)

GEDI.long<-melt(GEDI, id.vars="shot_number", measure.vars = colnames(GEDI)[-c(11,12,14, 16)])

ggplot(GEDI.long,
       aes(x=value))+
  geom_density(fill="red", color="black")+
  facet_wrap(~variable, scales="free")+
  theme_bw()

library(data.table)
GEDI.sub<-GEDI[GEDI$beam==unique(GEDI$beam)[1],]
GEDI.sub<-GEDI.sub[GEDI.sub$delta_time %between% list(quantile(GEDI.sub$delta_time,0.505),
                                              quantile(GEDI.sub$delta_time,0.507)),]


ggplot(GEDI.sub[order(GEDI.sub$delta_time),],
       aes(x=delta_time,y=elev_lowestmode))+
  
  geom_ribbon(aes(ymin=elev_lowestmode,ymax=rh100+elev_lowestmode),
              fill="forestgreen")+
  geom_path(aes(y=rh100+elev_lowestmode), color="forestgreen")+
  geom_path(color="black", size=0.1)+
  theme_bw()+theme(panel.grid = element_blank())

ggplot(GEDI.sub[order(GEDI.sub$delta_time),],
       aes(x=delta_time,y=rh100))+
  geom_ribbon(aes(ymin=0,ymax=rh100),
              fill="forestgreen", alpha=0.5)+
  geom_path(color="forestgreen")+
  theme_bw()+theme(panel.grid = element_blank())

library(leaflet)

pal <- colorNumeric(
  palette = "Blues",
  domain = GEDI.sub$rh100)

m1<-leaflet() %>%
  addCircleMarkers(GEDI.sub$lon,
                   GEDI.sub$lat,
                   radius = 1,
                   opacity = 1,
                   color = pal(GEDI.sub$rh100))  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)  %>%
  addLegend(colors = "blue", labels= c("GEDI Shots"),
            title ="GEDI Level2B") 
m1



library(rnaturalearth)
library(rGEDI)
library(rgdal)
library(raster)

country<-rnaturalearth::ne_countries(country="Bangladesh")
plot(country)
ext_bb<-extent(country)

# Study area boundary box coordinates
ul_lat<- ext_bb[c(4)]
lr_lat<- ext_bb[c(3)]
ul_lon<- ext_bb[c(1)]
lr_lon<- ext_bb[c(2)]

# Specifying the date range
daterange=c("2018-01-01","2020-06-24")

gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)

gedi_ls<-gLevel2B

gedi_ls<-as.vector(substr(gLevel2B, 
                          nchar("https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001/2019.10.02/")+1, 
                          200))


GEDI_files<-list.files("D:/GEDI/output", full.names = TRUE, pattern = ".txt")

gedi_ls_sub<-GEDI_files[(gsub("txt","h5",list.files("D:/GEDI/output", pattern="txt")) %in% gedi_ls)]

#What are the variables we can choose from?
vars<-colnames(data.table::fread(GEDI_files[1], nrow=1))[-(1:2)]
print(vars)

#choose the variables you want
output_variables<-vars[c(1,3)]

#Name the output file based on the project and variables included
Q1_output_file<-paste("D:/GEDI/output/rasters/Bangladesh_",
                      paste(output_variables, collapse = "_"),".txt", sep="")


library(sp)

for (i in 1:length(gedi_ls_sub)) {
  GEDI_df<-na.omit(data.table::fread(gedi_ls_sub[i], 
                             sep=" ", 
                             select = c('lon','lat',output_variables,'q1'), 
                             showProgress = FALSE))
  GEDI_sp<-SpatialPointsDataFrame(GEDI_df[,1:2],data=GEDI_df)
  GEDI_sp<-crop(GEDI_sp, countries)
  
  if(is.null(GEDI_sp)) next else {
    GEDI_df<-GEDI_sp@data
  
  data.table::fwrite(GEDI_df[GEDI_df$q1==1,1:(ncol(GEDI_df)-1)], 
                     file = Q1_output_file,
                     sep = " ",
                     na = "NA",
                     append=TRUE)
  }
  
  print(paste0((i/length(gedi_ls_sub))*100,' %'))
} 

GEDI<-fread(Q1_output_file)

GEDI_sp<-SpatialPointsDataFrame(coords=cbind(GEDI[,1],GEDI[,2]),
                                data=GEDI)

r<-raster(ext=extent(country), resolution=0.01)

library(raster)
GEDI_r<-rasterize(GEDI_sp, r, field="elev_lowestmode",fun=function(x, na.rm=TRUE) mean(x))
GEDI_rh100<-rasterize(GEDI_sp, r, field="rh100",fun=function(x, na.rm=TRUE) mean(x))


plot(GEDI_r, col=viridis::inferno(250))
plot(country, add=TRUE)

GEDI_rh100[GEDI_rh100>20]<-20
plot(GEDI_rh100,col=viridis::viridis(250))
plot(country, add=TRUE)







library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

option_ls<-data.frame(files=gedi_ls_sub,
     vars=t(output_variables),
     output=Q1_output_file, stringsAsFactors = FALSE)

parLapply(cl, option_ls[1:4,], function(x){
  vars<-t(x[,-c(1,ncol(x))])
  GEDI_df<-data.table::fread(x$files, 
                             sep=" ", 
                             select = c('lon','lat',vars,'q1'), 
                             showProgress = FALSE)
  data.table::fwrite(GEDI_df[GEDI_df$q1==1,1:(ncol(GEDI_df)-1)], 
                     file = x$output,
                     sep = " ",
                     na = "NA",
                     append=TRUE)
})

stopCluster(cl)
