library(raster)
library(rgeos)
library(rgdal)
library(gdalUtils)
library(gdata)
library(sf)

import_single_bgt <- function(aws.df = AWS.df, aws_name, sensor_name, addition = "", radius = 150, delete_raw_gmls = TRUE, redownload = FALSE){
  CleanGlobEnvir <- function(pattern){
    rm(list = ls(envir=globalenv())[
      grep(pattern, ls(envir=globalenv()))], envir = globalenv())
  }
  tryCatch({
    aws_name_trim <- getAWS_name_trim(aws.df = aws.df, aws_name = aws_name, addition = addition)
    if(file.exists(paste("data/BGT/", aws_name_trim, paste0("BGT_", aws_name_trim, ".shp"), sep="/")) & redownload == TRUE){  
      file.remove(paste("data/BGT/", aws_name_trim, paste0("BGT_", aws_name_trim, ".shp"), sep="/"))
      file.remove(paste("data/BGT/", aws_name_trim, paste0("BGT_", aws_name_trim, ".shx"), sep="/"))
      file.remove(paste("data/BGT/", aws_name_trim, paste0("BGT_", aws_name_trim, ".prj"), sep="/"))
      file.remove(paste("data/BGT/", aws_name_trim, paste0("BGT_", aws_name_trim, ".dbf"), sep="/"))
      warning("BGT shapefile for this aws is already found. They were removed and redownloaded.")
      
      import_single_bgt(aws.df = aws.df, aws_name = aws_name, sensor_name = sensor_name, radius = radius, delete_raw_gmls = delete_raw_gmls, redownload = TRUE)
    } else if(file.exists(paste("data/BGT/", aws_name_trim, paste0("BGT_", aws_name_trim, ".shp"), sep="/")) & redownload == FALSE){
      warning("BGT shapefile for this aws is already found and will used. Please remove it if you want to redownlolad it, or set redownload to TRUE..")
    } else {
      
      #aws station
      selected_aws <- select_single_aws(aws.df = aws.df, aws_name = aws_name, sensor_name = sensor_name) 
    
      #create 200 meter buffer arouond sensor
      surroundingBuffer <- buffer(selected_aws[["aws_rd.sp"]],width=radius)
      
      #get BBOX extent of buffer area
      surroundingExtent <- extent(surroundingBuffer)
      
      #create WFS string
      createWFS_string <- function(WFSbase_url, epsg, bbox){
        wfs_string <- paste("WFS:", WFSbase_url, "&SRSNAME=", epsg, "&BBOX=", bbox, sep="")
      return (wfs_string)}
      
      bgt_wfs <- createWFS_string(WFSbase_url = "https://geodata.nationaalgeoregister.nl/beta/bgt/wfs?request=getcapabilities" 
                      , epsg = "EPSG:28992"
                      , bbox = paste(toString(surroundingExtent@xmin), toString(surroundingExtent@ymin), toString(surroundingExtent@xmax), toString(surroundingExtent@ymax), sep=",")
                      )
      #xls_location <- paste(home_directory, wmoSC,data-raw, " bgt_names_features.xlsx", sep= "/")
      #bgt_objects.xlsx <- read.xls(xls_location[setup_settings], sheet="bgt_objects")
      
      bgt_objects.csv <- read.table(paste("wmoSC", "data-raw", "bgt_objects.csv", sep="/"), header = TRUE, sep=",",blank.lines.skip = TRUE)
      bgt_objects_shortname_list <- as.vector(data.frame(lapply(bgt_objects.csv[1], as.character), stringsAsFactors=FALSE)[,1])
      bgt_objects_name_list <- as.vector(data.frame(lapply(bgt_objects.csv[2], as.character), stringsAsFactors=FALSE)[,1])
      #bgt_features_perObject.xlsx <- read.xls(xls_location, sheet="bgt_features_perObject")
      bgt_features_perObject.csv <- read.table(paste("wmoSC", "data-raw", "bgt_features_perObject.csv", sep= "/"), header = TRUE, sep=",", blank.lines.skip = TRUE)
      
      bgt_begroeidterrein_features <- data.frame(lapply(bgt_features_perObject.csv[1], as.character), stringsAsFactors=FALSE)
      bgt_onbegroeidterrein_features <- data.frame(lapply(bgt_features_perObject.csv[2], as.character), stringsAsFactors=FALSE)
      bgt_waterdeel_features <- data.frame(lapply(bgt_features_perObject.csv[3], as.character), stringsAsFactors=FALSE)
      bgt_wegdeel_features <- data.frame(lapply(bgt_features_perObject.csv[4], as.character), stringsAsFactors=FALSE)
      bgt_scheiding_features <- data.frame(lapply(bgt_features_perObject.csv[5], as.character), stringsAsFactors=FALSE)
      bgt_pand_features <- data.frame(lapply(bgt_features_perObject.csv[6], as.character), stringsAsFactors=FALSE)
      bgt_overigbouwwerk_features <- data.frame(lapply(bgt_features_perObject.csv[7], as.character), stringsAsFactors=FALSE)
      bgt_spoor_features <- data.frame(lapply(bgt_features_perObject.csv[8], as.character), stringsAsFactors=FALSE)
      
      bgt_functioneelgebied_features <- data.frame(lapply(bgt_features_perObject.csv[9], as.character), stringsAsFactors=FALSE)
      
      #bgt_features_list <- list(bgt_begroeidterrein_features, bgt_onbegroeidterrein_features, bgt_waterdeel_features, bgt_wegdeel_features, bgt_scheiding_features, bgt_pand_features, bgt_overigbouwwerk_features, bgt_spoor_features, bgt_functioneelgebied_features)
      
      # omit <- function(features){
      #   #features[!apply(data == "", 1, all),]
      #   featureset <- na.omit(features)
      #   print(featureset)
      #   #features[complete.cases(features),]
      #   #return(features)
      # }
      # 
      # lapply(bgt_features_list, omit)
      # 
      # na.omit(bgt_waterdeel_features)
      bgt_begroeidterrein_features[complete.cases(bgt_begroeidterrein_features),]
      
      bgt_colNames <- c("bgt-functi", "bgt-fysiek", "bgt-type", "bronhouder", "gml_id", "naam", "object_type", "plus-fun_1", "plus-fun_2", "plus-funct", "plus-fys_2", "plus-fysie", "plus-sta_1", "plus-sta_2", "plus-statu", "plus-typ_1", "plus-type", "relatieveh", "status")
      bgt.df <- setNames(data.frame(matrix(ncol = 19, nrow = 0), stringsAsFactors = FALSE), bgt_colNames)
      bgt.df <- data.frame(lapply(bgt.df, function(x) if(is.logical(x)) {return(as.character(x))} else { return(x)}), stringsAsFactors=FALSE)
      
      #bgt_shape <- assign(paste("BGT", station, sep="_"), st_sfc(st_polygon(list())))
      
      #crs(BGT_deBilt) <- crs(epsg_rd)
      
      
      #adjust columns
      adjustColumns <- function(raw_data_location, object_name_short){
        print("adding columns...")
        shp <- readOGR(dsn = raw_data_location, layer = object_name_short, stringsAsFactors=FALSE)
        print(paste("before", length(names(shp))))
        ## add and remove the different columns to proceed later with correct rbind
        print(paste("in: ", object_name_short))
        if(object_name_short == "begroeidterreindeel"){
          shp$plus.funct <- NA
          shp$plus.fun_1 <- NA
          shp$plus.fun_2 <- NA
          shp$bgt.functi <- NA
          shp$plus.typ_1 <- NA
          shp$plus.type <- NA
          shp$bgt.type <- NA
          shp$naam <- NA
          shp$object_typ <- "begroeidterreindeel" 
          shp$object <- "vegetation"
          shp$colour <- "green"
          shp$hexColour <- "#5dcb1e"
          
          drops<- c('status_lee', 'status_cod', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'optalud_le', 'fysiekvoor', 'creationda', 'terminat_1', 'plus.fys_1', 'optalud', 'lv.publica', 'tijdstipre', 'kruinlijn_')
          shp <- shp[,!(names(shp) %in% drops)]
        
          } else if(object_name_short == "onbegroeidterreindeel"){
          shp$plus.funct <- NA
          shp$plus.fun_1 <- NA
          shp$plus.fun_2 <- NA
          shp$bgt.functi <- NA
          shp$plus.typ_1 <- NA
          shp$plus.type <- NA
          shp$bgt.type <- NA
          shp$naam <- NA
          shp$object_typ <- "onbegroeidterreindeel" 
          shp$object <- "barren"
          shp$colour <- "yellow"
          shp$hexColour<-"#f1f4c7"
          
          drops<- c('status_lee', 'status_cod', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'optalud_le', 'fysiekvoor', 'creationda', 'terminat_1', 'plus.fys_1', 'optalud', 'lv.publica', 'tijdstipre', 'kruinlijn_')
          shp <- shp[,!(names(shp) %in% drops)]
          
        } else if(object_name_short == "waterdeel"){
          shp$plus.fysie <- NA
          shp$plus.funct <- NA
          shp$plus.fun_1 <- NA
          shp$plus.fun_2 <- NA
          shp$bgt.functi <- NA
          shp$bgt.fysiek <- NA
          shp$plus.fys_2 <- NA
          shp$naam <- NA
          shp$object_typ <- "waterdeel" 
          shp$object <- "water"
          shp$colour <- "blue"
          shp$hexColour <- "#25e0ed"
          
          drops<- c('type_codes', 'status_lee', 'status_cod', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'creationda', 'terminat_1', 'lv.publica', 'tijdstipre', 'plus.type_')
          shp <- shp[,!(names(shp) %in% drops)]
          
        } else if(object_name_short == "wegdeel"){
          shp$plus.typ_1 <- NA
          shp$plus.type <- NA
          shp$bgt.type <- NA
          shp$naam <- NA
          shp$object_typ <- "wegdeel" 
          shp$object <- "road"
          shp$colour <- "gray"
          shp$hexColour <- "#93969b"
          
          drops<- c('status_lee', 'status_cod', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'plus.fys_1', 'lv.publica', 'kruinlijn_', 'tijdstipre', 'optalud', 'creationda', 'optalud_le', 'fysiekvoor', 'terminat_1', 'functie_co')
          shp <- shp[,!(names(shp) %in% drops)]
          
        } else if(object_name_short == "scheiding"){
          shp$plus.fysie <- NA
          shp$plus.funct <- NA
          shp$plus.fun_1 <- NA
          shp$plus.fun_2 <- NA
          shp$bgt.functi <- NA
          shp$bgt.fysiek <- NA
          shp$plus.fys_2 <- NA
          shp$naam <- NA
          shp$object_typ <- "scheiding" 
          shp$object <- "seperation"
          shp$colour <- "brown"
          shp$hexColour <- "#9b4b26"
          
          drops<- c('type_codes', 'status_lee', 'status_cod', 'plus.type_', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'creationda', 'terminat_1', 'lv.publica', 'tijdstipre')
          shp <- shp[,!(names(shp) %in% drops)]
          
        } else if(object_name_short == "pand"){
          shp$plus.fysie <- NA
          shp$plus.funct <- NA
          shp$plus.fun_1 <- NA
          shp$plus.fun_2 <- NA
          shp$bgt.functi <- NA
          shp$bgt.fysiek <- NA
          shp$plus.fys_2 <- NA
          shp$plus.typ_1 <- NA
          shp$plus.type <- NA
          shp$bgt.type <- NA
          shp$naam <- NA
          shp$object_typ <- "pand" 
          shp$object <- "building"
          shp$colour <- "red"
          shp$hexColour<- "#f0495f"
          
          drops<- c('identifica', 'status_lee', 'status_cod', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'creationda', 'terminat_1', 'lv.publica', 'tijdstipre')
          shp <- shp[,!(names(shp) %in% drops)]
          
        } else if(object_name_short == "overigbouwwerk"){
          shp$plus.fysie <- NA
          shp$plus.funct <- NA
          shp$plus.fun_1 <- NA
          shp$plus.fun_2 <- NA
          shp$bgt.functi <- NA
          shp$bgt.fysiek <- NA
          shp$plus.fys_2 <- NA
          shp$plus.typ_1 <- NA
          shp$plus.type <- NA
          shp$bgt.type <- NA
          shp$naam <- NA
          shp$object_typ <- "overigbouwwerk" 
          shp$object <- "other construction"
          shp$colour <- "black"
          shp$hexColour <- "#000000"
          
          drops<- c('type_codes', 'plus.type_', 'identifica', 'status_lee', 'status_cod', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'creationda', 'terminat_1', 'lv.publica', 'tijdstipre')
          shp <- shp[,!(names(shp) %in% drops)]
          
        } else if(object_name_short == "spoor"){
          shp$plus.fysie <- NA
          shp$bgt.fysiek <- NA
          shp$plus.fys_2 <- NA
          shp$plus.typ_1 <- NA
          shp$plus.type <- NA
          shp$bgt.type <- NA
          shp$naam <- NA
          shp$object_typ <- "spoor" 
          shp$object <- "railway"
          shp$colour <- "other"
          shp$hexColour <- "#ffffff"
          
          drops<- c('status_lee', 'status_cod', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'functie_co', 'creationda', 'terminat_1', 'lv.publica', 'tijdstipre')
          shp <- shp[,!(names(shp) %in% drops)]
          
        } else if(object_name_short == "functioneelgebied"){
          shp$plus.fysie <- NA
          shp$plus.funct <- NA
          shp$plus.fun_1 <- NA
          shp$plus.fun_2 <- NA
          shp$bgt.functi <- NA
          shp$bgt.fysiek <- NA
          shp$plus.fys_2 <- NA
          shp$plus.typ_1 <- NA
          shp$plus.type <- NA
          shp$object_typ <- "functioneelgebied" 
          shp$object <- "functional area"
          shp$colour <- "other"
          shp$hexColour <- "#ffffff"
          
          drops<- c('type_codes', 'status_lee', 'status_cod', 'plus.type_', 'inonderzoe', 'eindregist', 'terminatio', 'lokaalid', 'inonderz_1', 'creationda', 'terminat_1', 'naam_leeg', 'lv.publica', 'tijdstipre')
          shp <- shp[,!(names(shp) %in% drops)]
          
        }
        shp[ , order(names(shp))]
        print(paste("after", length(names(shp))))
        print(names(shp))
        writeOGR(obj = shp, dsn = raw_data_location, layer = object_name_short, driver = "ESRI Shapefile", overwrite_layer = TRUE)
        return (shp)
      }
      tmp.create_sp_atNext <<- FALSE
      read_bgt<-function(aws_name, wfs, bgt_object_name, object_name_short){
        bgt_directory <- paste("data", "BGT", sep="/")
        dir.create(paste(bgt_directory, aws_name_trim, sep="/"), showWarnings = FALSE)
        working_directory <- paste(bgt_directory, aws_name_trim, sep="/")
        dir.create(paste(working_directory, "raw", sep="/"), showWarnings = FALSE)
        #dir.create(paste(working_directory, "selections", sep="/"), showWarnings = FALSE)
        
        object_name <- paste("object", object_name_short, sep="_")
        object <- assign(object_name, SpatialPolygonsDataFrame(SpatialPolygons(list()), data=data.frame()))
        
        rowname <- paste(aws_name_trim,object_name_short,sep="_")
        print("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")
        print(rowname)
        
        raw_data_location <- paste(working_directory, "raw", sep="/")
        print(paste("raw data location: ", raw_data_location, sep=""))
        raw_shape_file_path <- paste(working_directory, "raw", "", sep="/")
        shape_file <- paste( raw_shape_file_path, object_name_short, ".shp", sep="")
        print(paste("shapefile: ", shape_file, sep=""))
    
        object <- tryCatch({
          ogr2ogr(src_datasource_name = wfs, dst_datasource_name = shape_file, layer = bgt_object_name, overwrite = TRUE)
          
          object_name <- readOGR(dsn = raw_data_location, layer = object_name_short, stringsAsFactors=FALSE)
          #object_name@data$object_type<- c(object_name_short)
          crs(object_name) <- crs(epsg_rd)
          if(class(object_name) == "SpatialLinesDataFrame"){
            lines.sf <- st_as_sf(object_name)
            polygons.sf <- st_polygonize(lines.sf)
            object_name <- sf:::as_Spatial(polygons.sf) 
          }
          #adjust columns
          object_name <- adjustColumns(raw_data_location, object_name_short)
          print(paste("field count object name: ", length(names(object_name@data))))
      
          crs(object_name) <- crs(epsg_rd)
          
          tmp.bgt_counter <<- tmp.bgt_counter + 1
          print(paste("counter pos. 1:", tmp.bgt_counter))
          
          if((tmp.bgt_counter == 1 & tmp.create_sp_atNext == FALSE) | (tmp.bgt_counter > 1 & tmp.create_sp_atNext == TRUE)){
            print("creating BGT_station.sp...")
            tmp.BGT_station.sp <<- object_name
            crs(tmp.BGT_station.sp) <- epsg_rd
            print(paste("length object:",length(object_name), sep=" "))
            print(paste("length BGT:",length(tmp.BGT_station.sp), sep=" "))
            tmp.create_sp_atNext <<- FALSE
            #print(length(BGT_shape))
          } else {  #print("in >1")
            print(paste("length new object:",length(object_name), sep=" "))
            print(paste("length BGT:",length(tmp.BGT_station.sp), sep=" "))
            
            tmp.BGT_station.sp <<- rbind(tmp.BGT_station.sp, object_name, makeUniqueIDs = TRUE)
            print(paste("new length BGT:",length(tmp.BGT_station.sp), sep=" "))
          }
          object_count <- length(object_name)
          
          entry_t <- data.frame(aws_name_trim, sensor_name, radius, object_name_short, TRUE, object_count, stringsAsFactors=FALSE)
          names(entry_t) <- c("AWS", "sensor_name", "radius", "object_type", "has_features", "feature_count")
          rownames(entry_t) <- rowname
          tmp.all_bgt_objects <<- rbind(tmp.all_bgt_objects, entry_t, stringsAsFactors=FALSE)
          
          count_entry <- data.frame(aws_name_trim, sensor_name, radius, object_name_short, object_count, stringsAsFactors = FALSE)
          names(count_entry) <- c("AWS", "sensor_name", "radius", "object_name", "object_count")
          tmp.objects_count_aws <<- rbind(tmp.objects_count_aws, count_entry, stringsAsFactors=FALSE)
        }, error=function(e){
          print(e)
          tmp.bgt_counter <<- tmp.bgt_counter + 1
          if(exists("tmp.BGT_station.sp") == FALSE){
            print("no exist")
            tmp.create_sp_atNext <<- TRUE
          }
          object <- NA
          no_count_entry <- data.frame(aws_name, sensor_name, radius, object_name_short, 0, stringsAsFactors = FALSE)
          names(no_count_entry) <- c("AWS", "sensor_name", "radius", "object_name", "object_count")
          print(no_count_entry)
          tmp.objects_count <<- rbind(tmp.objects_count_aws, no_count_entry, stringsAsFactors=FALSE)
          
          #remove empty shp
          shp <- paste(raw_shape_file_path, object_name_short, ".shp", sep="")
          dbf <- paste(raw_shape_file_path, object_name_short, ".dbf", sep="")
          prj <- paste(raw_shape_file_path, object_name_short, ".prj", sep="")
          shx <- paste(raw_shape_file_path, object_name_short, ".shx", sep="")
          #file.remove(shp, dbf, prj, shx)
          
          print(e)
          #enter that object has no features
          entry_f <- data.frame(aws_name_trim, sensor_name, radius, object_name_short, FALSE, 0, stringsAsFactors=FALSE)
          names(entry_f) <- c("AWS", "sensor_name", "radius", "object_type", "has_features", "feature_count")
          rownames(entry_f) <- rowname
          tmp.all_bgt_objects <<- rbind(tmp.all_bgt_objects, entry_f, stringsAsFactors=FALSE)
          message(paste("No features found for", object_name_short, "in", aws_name, sep=" "))
          object_name <- NA
        }
        )
        
        print(paste("counter pos. 2:" ,tmp.bgt_counter))
        print(paste("total object types count", length(bgt_objects_name_list)))
        if(tmp.bgt_counter == length(bgt_objects_name_list)){
          shp_name <- paste("BGT", aws_name_trim, sep="_")
          tmp.BGT_station.sp$AREA <-sapply(slot(tmp.BGT_station.sp, 'polygons'), function(i) slot(i, 'area')) 
          proj4string(tmp.BGT_station.sp) <- epsg_rd
          writeOGR(obj = tmp.BGT_station.sp, dsn = working_directory, layer = shp_name, driver = "ESRI Shapefile", overwrite_layer = TRUE)
          tmp.BGT_station.sf <<- st_as_sf(tmp.BGT_station.sp)
          st_crs(tmp.BGT_station.sf, epsg_rd)
          tmp.bgt_counter <<- 0
          return(tmp.BGT_station.sf)
        }
      }
      tmp.bgt_counter <<- 0
      #read BGT and create one BGT shp for the aws.
      tmp.objects_count_aws <<- data.frame(aws_name = as.character(), sensor_name = as.character(), radius = as.numeric(), objectname = as.character(), count = as.numeric(), stringsAsFactors=FALSE)
      objects_counts_names <- c("AWS", "sensor_name", "radius", "object_name", "object_count")
      names(tmp.objects_count_aws) <- objects_counts_names
      
      tmp.all_bgt_objects <<- data.frame(AWS = character(0), sensor_name = character(0), radius = numeric(0), feature_type = character(0), has_features = logical(0), features_count = numeric(0), stringsAsFactors = FALSE) 
     
      mapply(read_bgt,aws_name = aws_name, wfs = bgt_wfs, bgt_object_name = bgt_objects_name_list, object_name_short = bgt_objects_shortname_list)
      
      BGT_aws.sf <- tmp.BGT_station.sf
      st_crs(BGT_aws.sf, epsg_rd)
      objects_count <- tmp.objects_count_aws 
      all_bgt_objects <- tmp.all_bgt_objects 
      CleanGlobEnvir(pattern = "tmp")
      if(delete_raw_gmls == TRUE){
       unlink(paste("data", "BGT", aws_name_trim, "raw" ,sep="/"), recursive = TRUE)
      }
      return(list("BGT.sf" = BGT_aws.sf, "relevant_objects_count"=objects_count, "all_objects_count"=all_bgt_objects))
    }
  }, error=function(e){
    print(e)
    CleanGlobEnvir(pattern = "tmp")
  })
}
