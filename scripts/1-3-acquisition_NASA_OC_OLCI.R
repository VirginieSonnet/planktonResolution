###########################
# NetCDF files processing #
###########################

# Date: 2022-12-12
# Author: Virginie Sonnet 

# Goal: Process all the netcdf files from OLCI downloaded from the NASA Ocean Color website 
# for Narragansett Bay.

# Outputs: 
# - raw_gsoOptics_OLCI_metadata.csv
# - raw_gsoOptics_OLCI_chla_3x4_subset.csv

# make sure that no object is in the environment and collect the garbage 
rm(list=ls())
gc()

library(tidyverse)
library(lubridate)
library(ncdf4)


# Files -------------------------------------------------------------------

# get all the filenames 
files = tibble(filepath=list.files(path="data/olci_oc",pattern="*.nc",
                                   full.names = TRUE))


# extract the dates 
files <- files %>% 
    mutate(filename=str_split_i(filepath,"/",3),
           sampling_date=str_split_i(filename,"[.]",2),
           day=str_split_i(sampling_date,"T",1),
           time = as.numeric(str_split_i(sampling_date,"T",2)),
           satellite=str_split_i(filename,"_",1)) %>% 
    select(-sampling_date) %>% 
    # create the UTC date 
    mutate(UTC=ymd_hms(paste(day,time,sep=""))) %>% 
    # convert the UTC date to EST_EDT
    mutate(EST_EDT=with_tz(UTC, tzone="America/New_York"))


# create a set of unique ids: swath_id and sampling_id (group swaths from the
# same satellite and same day)
files <- files %>% 
    # sort by date 
    arrange(UTC) %>% 
    # add the unique swath_id
    mutate(swath_id=row_number()) %>% 
    # add the sampling_id based on satellite and day
    group_by(satellite,day) %>% 
    mutate(sampling_id=cur_group_id()) %>% 
    ungroup()



# Loop over files ---------------------------------------------------------

## Create 2 empty tibbles for storage
# metadata
data <- c()
# selected pixels with chlorophyll a and latitude
selected <- c()

for (i in files$swath_id){
    
    metadata <- filter(files,swath_id==i)
    print(paste("Processing:",metadata$filepath))
    
    ## Extract the variables 
    # open the nc file connection 
    nc <- nc_open(metadata$filepath)
    
    # get the variables 
    flags <- ncvar_get(nc, "geophysical_data/l2_flags") %>% as.vector()
    chl <- ncvar_get(nc, "geophysical_data/chlor_a") 
    lat <- ncvar_get(nc, "navigation_data/latitude") 
    lon <- ncvar_get(nc, "navigation_data/longitude") 
    
    # extract the dimensions 
    nrows = nc$dim$pixels_per_line$len
    ncols = nc$dim$number_of_lines$len
    
    # close the nc file connection
    nc_close(nc)
    
    # create a tibble
    swath <- tibble(chl=as.vector(chl),
                    lat=as.vector(lat),
                    lon=as.vector(lon),
                    flags) %>% 
        # add a variable to define land pixels
        mutate(land=ifelse(flags %in% c(1073741826,1073742082,1073742338),TRUE,FALSE)) %>% 
        # add an id for each pixel 
        mutate(pixel_id=row_number())
    
    
    ## Add to the metadata 
    # average chlorophyll a 
    metadata$avg_chl <- mean(chl,na.rm=TRUE) 
    
    # percentage of unresolved 
    notland <- filter(swath,land==FALSE)
    tot <- nrow(notland)
    metadata$p_unresolved <- nrow(filter(notland,is.na(chl)))/tot*100
    
    
    ## Extract the 3x4 pixels 
    # calculate how far are each pixels from the GSO Dock
    square = swath %>% add_column(GSOlat=41.492235, GSOlon=-71.418863) %>% 
        mutate(difflat = (GSOlat-lat)^2,
               difflon=(GSOlon-lon)^2) %>% 
        mutate(dist=sqrt(difflat+difflon)) %>% 
        arrange(dist)
    
    metadata$closest <- square$dist[1]
    if (square$dist[1] > 0.005){
        print("GSO Dock is not included in this swath.")
        metadata$GSO <- FALSE 
        
    } else {
        
        metadata$GSO <- TRUE
        
        # selected closest pixel 
        latgoal=square$lat[1]
        longoal=square$lon[1]
        
        # pixel_id corresponding to the closest pixel 
        match=filter(swath,lat==latgoal,lon==longoal)
        idx=max(match$pixel_id) # take the max just in case 2 pixels have same coordinates
        
        # calculate the indices corresponding to the closest pixel 
        colid = ceiling(idx/nrows)
        rowid = idx-((colid-1)*nrows)
        
        # check that the indices found do correspond to our targeted pixel
        if (!(latgoal == lat[rowid,colid] & longoal == lon[rowid,colid])){
            print("There is a problem with the latitude/longitude matching.")
            metadata$matching <- "wrong"
            
        } else {
            
            metadata$matching <- "correct"
            
            # 3x4 square east of it 
            sq_row <- rep(c(rowid+1,rowid+2,rowid+3,rowid+4),each=3)
            sq_col <- rep(c(colid-1,colid,colid+1),times=4)
            pixel_id=(sq_col-1)*nrows+sq_row
            ids=c(idx,pixel_id)
            
            # select the data
            swath <- filter(swath,pixel_id %in% ids) %>% 
                mutate(swath_id=i,
                       type=ifelse(lat==latgoal & lon==longoal,"GSO","NBay")) 
            
            # add to the previous data 
            selected <- bind_rows(selected,swath)
        }

    }
    data <- bind_rows(metadata,data)
}



# Export ------------------------------------------------------------------

# selected chlorophyll data 
write_csv(data, "data/raw_gsoOptics_OLCI_metadata.csv")
write_csv(selected,"data/raw_gsoOptics_OLCI_chla_3x4_subset.csv")
