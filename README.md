---
title: "Introduction to GEDI Processing"
output:
  html_document: default
  pdf_document: default
---

# Getting Started with Finding Available GEDI (from rGEDI tutorial page)

## Installation of rGEDI
```r
#The CRAN version:
install.packages("rGEDI")

#The development version:
library(devtools)
devtools::install_git("https://github.com/carlos-alberto-silva/rGEDI", dependencies = TRUE)

# loading rGEDI package
library(rGEDI)

```    

## Find GEDI data within your study area (GEDI finder tool)

```r
# Study area boundary box coordinates
ul_lat<- -44.0654
lr_lat<- -44.17246
ul_lon<- -13.76913
lr_lon<- -13.67646
```
Better yet, I like to import a raster or shapefile of the area of interest to grab the bounding coordinates like this:

```r
#pick your country of interest
countries<-rnaturalearth::ne_countries()

#read in a shapefile of interest
country<-readOGR("D:/US_conterminous_extent.shp")

#extent to get GEDI data
ext_bb<-extent(country)

# Study area boundary box coordinates (less manual. Woo!)
ul_lat<- ext_bb[c(4)]
lr_lat<- ext_bb[c(3)]
ul_lon<- ext_bb[c(1)]
lr_lon<- ext_bb[c(2)]
```

Note the date range of the data you desire. This might align with growing season or some other event of interest.
```r
# Specifying the date range
daterange=c("2019-07-01","2020-05-22")
```

Find available GEDI data:
```r
# Get path to GEDI data
gLevel1B<-gedifinder(product="GEDI01_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2A<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
```
## Downloading GEDI data

Using rGEDI (this is slow!)
```r
# Set output dir for downloading the files
outdir=getwd()

# Downloading GEDI data
gediDownload(filepath=gLevel1B,outdir=outdir)
gediDownload(filepath=gLevel2A,outdir=outdir)
gediDownload(filepath=gLevel2B,outdir=outdir)
```
The method I have adopted uses Cygwin64 terminal to run a modified GEDI download bash script:
```r

gedi_ls<-gLevel2B

#this defines the ftp path, so we can sub in our local directory path
gedi_path<-as.vector(substr(gedi_ls,1,
                            nchar("https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001/2019.04.18/")))

#we modify the list so we know what we have already downloaded
gedi_ls<-as.vector(substr(gedi_ls, 
                          nchar("https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001/2019.04.18/")+1, 
                          200))

#remove GEDI files that are already downloaded to the local drive
gedi_ls_sub<-gedi_ls[!(gedi_ls %in% list.files("D:/GEDI/", pattern="h5"))]

#add the correct file path back to the files for downloading
gedi_ls_sub<-paste0(gedi_path[!(gedi_ls %in% list.files("D:/GEDI/", pattern="h5"))],gedi_ls_sub)

#write this file list
writeLines(gedi_ls_sub, "D:/GEDI/GEDI_ls_sub.txt")

```

Now (yes this is silly...) I paste the massive file list into the download bash script. Here's where using Nathan's method would be a winner if you could run it easily from R. I have tried automating the creation of this script without success. Maybe we can integrate these two methods for speed!


# Getting started preprocessing GEDI data

## Setup file

Now we will modify a simple setup file so we can blaze thorugh a bunch of GEDI data:

```{r}
library(h5)
source("R/GEDI_processing_FUN.R")
beams<-c("0000","0001","0010","0011",
         "0101","0110", "1000", "1011")

GEDI.dir<-"D:/GEDI/"
output.dir<-"D:/GEDI/output/"

file_ls<-subset.files(GEDI.dir, output.dir)
```

Here we want to modify the directory where we find the GEDI data (`GEDI.dir`) and where the new text files will be written (`output.dir`). `subset.files()` simply makes sure you are not re-processing h5 GEDI files that already exist.

## GEDI data structure
First, let's take a look 
```{r echo=TRUE}
print(file_ls[1:4])
```

For an example of the GEDI data structure I'm going to choose a single orbit that overlaps with Gabon and Pongara National Park:

```{r echo=TRUE}
file_ls<-"E:/Gabon/GEDI/current_release_correct/GEDI02_B_2019111040155_O02008_T04616_02_001_01.h5"
```

Normally, this portion of the script is included in a loop or similar, but this is the basic method of reading h5 GEDI data in:

```{r}
i=1

file <- h5::h5file(name = file_ls[i],
                 mode = "a")

```

Now that we have a GEDI file in R we need to look at the stucture. It has a very complex hierarchical structure:

```{r}
  data.sets<-h5::list.datasets(file)

data.sets[1:180]

```


**WOW!** That is a lot of info to grab (180 categories X 8 beams!). So you need to choose what is important to you for processing.

Now, within each beam we can subset by the specific data we need. 

# GEDI Preprocessing
## Preprocessing Function
Ideally we want to do all of this automatically for lots of GEDI data, so we can use this simple function:

```r
preprocess.GEDI(input.file=x, 
                      filename=filename)
```

`preprocess.GEDI` takes h5 files, selects useful metrics, and outputs a simple txt file for further processing. This greatly simplifies the pipeline. We can modify this function to include specific information from the h5 file, but for now it includes: 

"lon"
"lat"
"elev_lowestmode"
"elev_TDX"
"rh100"
"cover"
"pai"
"fhd"
"num_detectedmodes"
"modis_treecover"
"q1" "surface_flag"
"sensitivity" 
"shot_number" "delta_time"
"beam"

I find these variables give me everything I need from the L2B data product.

## Process in parallel

Now we can use the `parallel` package to process everything with lightning speeds!

```r
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

parLapply(cl, file_ls, function(x){
  
  source("R/GEDI_setup.R")
  
  filename<-gsub(".h5",".txt", 
                 gsub(GEDI.dir,output.dir, x) )
  
  try(preprocess.GEDI(input.file=x, 
                      filename=filename), silent=TRUE)
})

stopCluster(cl)

```
Check the output directory to see the GEDI files being written to your drive.

## Read a processed GEDI .txt file

Now we can easily use the preprocessed GEDI data to do some science.

```{r echo=TRUE,  warning=FALSE}
#Get our processed GEDI files
GEDI.processed.files<-list.files(output.dir, pattern="txt", full.names = TRUE)

library(data.table)
#read a single orbit file in
GEDI<-data.table::fread(GEDI.processed.files[1])

#subset to "quality 1" shots (the good ones!)
GEDI<-GEDI[GEDI$q1==1,]

#load some some packages for viz
library(ggplot2)
library(reshape)

GEDI.long<-melt(GEDI, id.vars="shot_number", measure.vars = colnames(GEDI)[-c(11,12,14, 16)])

ggplot(GEDI.long,
       aes(x=value))+
  geom_density(fill="red", color="black")+
  facet_wrap(~variable, scales="free")+
  theme_bw()

```

With this we can see the distribution of values for each variable. Now we can visualize some GEDI transects:

```{r echo=FALSE}
#lets grab a small subset based on delta_time
GEDI.sub<-GEDI[GEDI$beam==unique(GEDI$beam)[1],]
GEDI.sub<-GEDI.sub[GEDI.sub$delta_time %between% list(quantile(GEDI.sub$delta_time,0.505),
                                              quantile(GEDI.sub$delta_time,0.507)),]

```

Let's take a look at the profiles with elevation and canopy height:

```{r echo=TRUE}

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
```

Where are we??? Let's look in a map:
```{r echo=FALSE}
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
```

# Gridding GEDI data

Now that we have a bunch of GEDI data processed we can get to work on gridding some GEDI variables. Let's focus on elevation and canopy height (rh100) for this tutorial.

First we decide where we want to grid the GEDI data...

In honor of the acceped paper, let's check out Bangladesh!

```{r echo=FALSE, warning=FALSE}

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

```

Here we just got a list of GEDI files available in Bangladesh and subset based on the files we have available locally.

Now, want to put all of the GEDI shots within Bangladesh in one file before we grid. First, we decide on variables and output file name:

```{r}

#What are the variables we can choose from?
vars<-colnames(data.table::fread(GEDI_files[1], nrow=1))[-(1:2)]
print(vars)

#choose the variables you want
output_variables<-vars[c(1,3)]

#Name the output file based on the project and variables included
Q1_output_file<-paste("D:/GEDI/output/rasters/Bangladesh_",
                      paste(output_variables, collapse = "_"),".txt", sep="")


```


Now, we iteratively open GEDI txt files, clip the orbits, and save as a single new file:

```r


library(sp)

for (i in 1:length(gedi_ls_sub)) {
  GEDI_df<-na.omit(data.table::fread(gedi_ls_sub[i], 
                             sep=" ", 
                             select = c('lon','lat',output_variables,'q1'), 
                             showProgress = FALSE))
  GEDI_sp<-SpatialPointsDataFrame(GEDI_df[,1:2],data=GEDI_df)
  GEDI_sp<-crop(GEDI_sp, country)
  
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


```

Now we can read our newly created file and take a look:

```{r echo=TRUE}

GEDI<-fread(Q1_output_file)

hist(GEDI$elev_lowestmode, breaks=100)
hist(GEDI$rh100, breaks=100)

```


Finally, we can easily convert the GEDI shots into spatial points and rasterize them to a grid resolution of our choice:

```{r echo=TRUE}

GEDI_sp<-SpatialPointsDataFrame(coords=cbind(GEDI[,1],GEDI[,2]),
                                data=GEDI)


library(raster)
r<-raster(ext=extent(country), resolution=0.01)
GEDI_r<-rasterize(GEDI_sp, r, field="elev_lowestmode",fun=function(x, na.rm=TRUE) mean(x))
GEDI_rh100<-rasterize(GEDI_sp, r, field="rh100",fun=function(x, na.rm=TRUE) mean(x))


plot(GEDI_r, col=viridis::inferno(250))
plot(country, add=TRUE)

GEDI_rh100[GEDI_rh100>20]<-20
plot(GEDI_rh100,col=viridis::viridis(250))
plot(country, add=TRUE)

```

