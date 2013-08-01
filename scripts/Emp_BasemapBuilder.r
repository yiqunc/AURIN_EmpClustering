# required packages
library(maptools)   # for geospatial services; also loads foreign and sp
library(rgdal)      # for map projection work; also loads sp
library(rgeos)
library(gdata)
library(RJSONIO)
library(PBSmapping) # for GIS_like geospatial object manipulation / anslysis including poly

# set up the working directory
setwd("/Users/Shared/Documents/AURIN/R")

div1 = read.table("./data/ABS_data_by_DZN/DZN/div1_v3.csv", header=TRUE, na.strings="NA", sep="\t", blank.lines.skip=TRUE, fill=TRUE, colClasses="character")
div2 = read.table("./data/ABS_data_by_DZN/DZN/div2_v3.csv", header=TRUE, na.strings="NA", sep="\t", blank.lines.skip=TRUE, fill=TRUE, colClasses="character")
div3 = read.table("./data/ABS_data_by_DZN/DZN/div3_v3.csv", header=TRUE, na.strings="NA", sep="\t", blank.lines.skip=TRUE, fill=TRUE, colClasses="character")
div4 = read.table("./data/ABS_data_by_DZN/DZN/div4_v3.csv", header=TRUE, na.strings="NA", sep="\t", blank.lines.skip=TRUE, fill=TRUE, colClasses="character")

# merge DZN polygons with Employment data via DZN codes
# collapse child jobcat numbers to parent jobcat in the original job data before splitting job numbers to the split-polygons
# This will cause a problem. 
# let's suppose that A000 contains 5 direct children categories: X0100, X0200, X0300, X0400, X0500
# and the zonecodes in A000 contains IN1Z and FZ, and the X0400 only contains FZ, and the rest children only contains IN1Z.
# In a DZN zone area, it happens to exist two split-polygons need to be assigned with job numbers. Split-polygon P1 belongs to IN1Z and its area is 90(m*m)
# Split-polygon P2 belongs to FZ, and its area is 10(m*m).
# Before collapse, suppose the job numbers in that DZN area are:
# A0000, X0100, X0200, X0300, x0400, X0500
#     2,    10,     0,     5,     5,     0
# 
# if we collapse in the original data, it turns to be:
# A0000, X0100, X0200, X0300, x0400, X0500
#    22,    10,     0,     5,     5,     0
# then we assign jobnumbers in each category on area proportion to a split-polygons regarding to its zonecode:
# for P1, which belongs to IN1Z, and for its jobnumbers in A000 jobcat will be 90%*22=19.8
# for P1, which belongs to IN1Z, and for its jobnumbers in X0100/X0200/X0300/x0400/X0500 will be 100%*10, 100%*0, 100%*5, 0%*5, 100%*0
# for P2, which belongs to FZ, and for its jobnumbers in A000 jobcat will be 10%*22=2.2
# for P2, which belongs to FZ, and for its jobnumbers in X0100/X0200/X0300/X0500 will be 0%*10, 0%*0, 0%*5, 100%*5, 0%*0
#     A0000, X0100, X0200, X0300, x0400, X0500
# P1   19.8,    10,     0,     5,     0,     0
# P2    2.2,     0,     0,     0,     5,     0
# So if we do it this way, in the split polygons, the sum of job number of children jobcats does not equal to the parent. 
# This is because the zonecode structure in children and parent jobcats are not the same.
# 
# Conclusion: we can build collapse jobnumbers for the original job data, but it cannot be used as a base to assign job numbers to split-polygons.
# To assign correct job numbers to split-polygons: (1) split job numbers directly on jobcat, (2) collapse in the end. e.g:
# A0000, X0100, X0200, X0300, x0400, X0500
#     2,    10,     0,     5,     5,     0
# for P1, which belongs to IN1Z, and for its jobnumbers in A000 jobcat will be 90%*2=1.8
# for P1, which belongs to IN1Z, and for its jobnumbers in X0100/X0200/X0300/x0400/X0500 will be 100%*10, 100%*0, 100%*5, 0%*5, 100%*0
# for P2, which belongs to FZ, and for its jobnumbers in A000 jobcat will be 10%*2=0.2
# for P2, which belongs to FZ, and for its jobnumbers in X0100/X0200/X0300/X0500 will be 0%*10, 0%*0, 0%*5, 100%*5, 0%*0
#     A0000, X0100, X0200, X0300, x0400, X0500
# P1    1.8,    10,     0,     5,     0,     0
# P2    0.2,     0,     0,     0,     5,     0
# then collapse:
#     A0000, X0100, X0200, X0300, x0400, X0500
# P1   16.8,    10,     0,     5,     0,     0
# P2    5.2,     0,     0,     0,     5,     0
# v1.2 collapses jobnumber before split. v1.3 splits jobnumber before collapse

f_datamerge <- function(){
  tmpdf = read.table("./data/ABS_data_by_DZN/DZN/AllCategories.txt", header=TRUE, na.strings="-", blank.lines.skip = TRUE, fill = TRUE)
  # 2800 is somehow missing in the AllCategories.txt data but exists in the job categories. Have to append this column for integrity
  tmpdf[,"X2800"] = 0
  tmpmx = as.matrix(tmpdf)
  tmpmx[which(is.na(tmpmx))] = 0
  
  empDataFrame = as.data.frame(tmpmx)
  
  # join attribute data to polygon data on dzn code, all polygon data rows are kept
  empPolyRawDataFrame <- readOGR(dsn="./outputs",layer="DZN_Original",encoding="utf8")
   
  empPolyRawDataFrame@data[,"orgSeqId"] = c(1:nrow(empPolyRawDataFrame@data))
  joinedDF = merge(empPolyRawDataFrame@data, empDataFrame, by.x = "Vicdznp06", by.y = "DZN", all.x=TRUE, sort=FALSE)
  
  orderedJoinedDF = joinedDF[order(joinedDF[,"orgSeqId"]), ]
  # replace all unmatched rows' job number with 0
  filter = is.na(orderedJoinedDF[,"FULL"])
  orderedJoinedDF[filter,c(13:length(colnames(orderedJoinedDF)))] = 0 
  empPolyRawDataFrame@data = orderedJoinedDF
  
  writeOGR(obj=empPolyRawDataFrame, dsn="./outputs", layer="DZN_X_Employment_Orginal", driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
}

f_datamerge_2011 <- function(){
  tmpdf = read.table("./data/ABS_data_by_DZN_2011/INDP_X_DZN_2011_V1.txt", sep=",", header=TRUE, na.strings="-", blank.lines.skip = TRUE, fill = TRUE)
  # 2800 is somehow missing in the AllCategories.txt data but exists in the job categories. Have to append this column for integrity
  tmpdf[,"X2800"] = 0
  tmpmx = as.matrix(tmpdf)
  tmpmx[which(is.na(tmpmx))] = 0
  
  empDataFrame = as.data.frame(tmpmx)
  
  # join attribute data to polygon data on dzn code, all polygon data rows are kept
  empPolyRawDataFrame <- readOGR(dsn="./data/ABS_data_by_DZN_2011",layer="DZN_2011_VIC",encoding="utf8")
  
  empPolyRawDataFrame@data[,"orgSeqId"] = c(1:nrow(empPolyRawDataFrame@data))
  joinedDF = merge(empPolyRawDataFrame@data, empDataFrame, by.x = "DZN_CODE11", by.y = "DZN", all.x=TRUE, sort=FALSE)
  
  orderedJoinedDF = joinedDF[order(joinedDF[,"orgSeqId"]), ]
  # replace all unmatched rows' job number with 0
  filter = is.na(orderedJoinedDF[,"DZN_CODE11"])
  orderedJoinedDF[filter,c(13:length(colnames(orderedJoinedDF)))] = 0 
  empPolyRawDataFrame@data = orderedJoinedDF
  
  writeOGR(obj=empPolyRawDataFrame, dsn="./data/ABS_data_by_DZN_2011", layer="DZN_X_Employment_2011_Orginal", driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
}

# based on Sophie's work, parsing the summarized "jobdivision-zonecode" table into 4 tables, which contains the "jobdivision-zonecode" relationship on each of 4 jobdivision levels respectively. 
f_parse_divzonecode <- function(){
  # parsing division-zonecode data level1
  divDataFrame = read.table("./data/ABS_data_by_DZN/DZN/divdata_v3.txt", header=TRUE, na.strings="", sep="\t", blank.lines.skip = TRUE, fill = TRUE)
  div1Filter = !is.na(divDataFrame[,"Div1"])
  div1DF = divDataFrame[div1Filter,c("Div1","Div2")]
  attr(div1DF, "names")[1] = "Div1Code"
  attr(div1DF, "names")[2] = "Div1Name"
  
  # parsing division-zonecode data level2
  div2 = divDataFrame[,c("Div1","Div2","Div3")]
  div2DF <- data.frame(Div1Code=NULL, Div2Code=NULL, Div2Name=NULL)
  div1Code = ""
  for(i in 1:nrow(div2)){
      if (!is.na(div2[i,"Div1"])){
      div1Code = as.character(div2[i,"Div1"])
      #print(sprintf("div1: %s", div1Code))
      next
    }
    
    if (!is.na(div2[i,"Div2"])){
      #print(sprintf("div2: %s", as.character(div2[i,"Div2"])))
      div2DF <- rbind(div2DF, data.frame(Div1Code=div1Code, Div2Code=as.character(div2[i,"Div2"]), Div2Name=as.character(div2[i,"Div3"])))
    }
}

  # parsing division-zonecode data level3
  div3 = divDataFrame[,c("Div1","Div2","Div3", "Div4")]
  div3DF <- data.frame(Div1Code=NULL, Div2Code=NULL, Div3Code=NULL, Div3Name=NULL)
  div1Code = ""
  div2Code = ""
  for(i in 1:nrow(div3)){
    if (!is.na(div3[i,"Div1"])){
    div1Code = as.character(div3[i,"Div1"])
    print(sprintf("div1: %s", div1Code))
    next
  }
  
  if (!is.na(div3[i,"Div2"])){
    div2Code = as.character(div3[i,"Div2"])
    print(sprintf("--div2: %s", div2Code))
    next
  }
  
  if (!is.na(div3[i,"Div3"]) & !is.na(as.numeric(as.character(div3[i,"Div3"])))){
    print(sprintf("----div3: %s", as.character(div3[i,"Div3"])))
    div3DF <- rbind(div3DF, data.frame(Div1Code=div1Code, Div2Code = div2Code, Div3Code=as.character(div3[i,"Div3"]), Div3Name=as.character(div3[i,"Div4"])))
  }
}

  # parsing division-zonecode data level4
  div4 = divDataFrame[,c("Div1","Div2","Div3", "Div4", "Div4Name", "PlanZoneCode")]
  div4DF <- data.frame(Div1Code=NULL, Div2Code=NULL, Div3Code=NULL, Div4Code=NULL, Div4Name=NULL, ZoneCodes=NULL)
  div1Code = ""
  div2Code = ""
  div3Code = ""
  for(i in 1:nrow(div4)){
    if (!is.na(div4[i,"Div1"])){
      div1Code = as.character(div4[i,"Div1"])
      print(sprintf("div1: %s", div1Code))
      next
    }
    
    if (!is.na(div4[i,"Div2"])){
      div2Code = as.character(div4[i,"Div2"])
      print(sprintf("--div2: %s", div2Code))
      next
    }
    
    if (!is.na(div4[i,"Div3"]) & !is.na(as.numeric(as.character(div4[i,"Div3"])))){
      div3Code = as.character(div4[i,"Div3"])
      print(sprintf("----div3: %s", div3Code))
      next
    }
    
    if (!is.na(div4[i,"Div4"]) & !is.na(as.numeric(as.character(div4[i,"Div4"])))){
      print(sprintf("------div4: %s", as.character(div4[i,"Div4"])))
      div4DF <- rbind(div4DF, data.frame(Div1Code=div1Code, Div2Code = div2Code, Div3Code = div3Code, Div4Code=as.character(div4[i,"Div4"]), Div4Name=as.character(div4[i,"Div4Name"]), ZoneCodes=div4[i,"PlanZoneCode"]))
    }
  }
  write.table(div4DF, file="./data/ABS_data_by_DZN/DZN/div4_v3.csv", row.names = FALSE, sep="\t")

  # build zonecode string for div3-div2-div1
  div3DF[,"ZoneCodes"] = ""
  for(k in 1: nrow(div3DF)){
    div3Code = as.character(div3DF[k,"Div3Code"])
    filter = div4DF[,"Div3Code"] == div3Code
    tmpDF = div4DF[filter,]
    zonecodes = ""
    for(i in 1:nrow(tmpDF)){
      if(!is.na(as.character(tmpDF[i, "ZoneCodes"]))){
        zonecodes = paste(zonecodes,as.character(tmpDF[i, "ZoneCodes"]), sep=",")
      }
    }
    if (zonecodes == ""){
      div3DF[k,"ZoneCodes"] = NA
    } else {
      vec = trim(strsplit(zonecodes,",")[[1]])
      vec = vec[!duplicated(vec) & vec!=""]
      div3DF[k,"ZoneCodes"] = paste(vec,collapse=",")
    }
  }
  write.table(div3DF, file="./data/ABS_data_by_DZN/DZN/div3_v3.csv", row.names = FALSE, sep="\t")
  
  div2DF[,"ZoneCodes"] = ""
  for(k in 1: nrow(div2DF)){
    div2Code = as.character(div2DF[k,"Div2Code"])
    filter = div4DF[,"Div2Code"] == div2Code
    tmpDF = div4DF[filter,]
    zonecodes = ""
    for(i in 1:nrow(tmpDF)){
      if(!is.na(as.character(tmpDF[i, "ZoneCodes"]))){
        zonecodes = paste(zonecodes,as.character(tmpDF[i, "ZoneCodes"]), sep=",")
      }
    }
    if (zonecodes == ""){
      div2DF[k,"ZoneCodes"] = NA
    } else {
      vec = trim(strsplit(zonecodes,",")[[1]])
      vec = vec[!duplicated(vec) & vec!=""]
      div2DF[k,"ZoneCodes"] = paste(vec,collapse=",")
    }
  }
  write.table(div2DF, file="./data/ABS_data_by_DZN/DZN/div2_v3.csv", row.names = FALSE, sep="\t")

  div1DF[,"ZoneCodes"] = ""
  for(k in 1: nrow(div1DF)){
    div1Code = as.character(div1DF[k,"Div1Code"])
    filter = div4DF[,"Div1Code"] == div1Code
    tmpDF = div4DF[filter,]
    zonecodes = ""
    for(i in 1:nrow(tmpDF)){
      if(!is.na(as.character(tmpDF[i, "ZoneCodes"]))){
        zonecodes = paste(zonecodes,as.character(tmpDF[i, "ZoneCodes"]), sep=",")
      }
    }
    if (zonecodes == ""){
      div1DF[k,"ZoneCodes"] = NA
    } else {
      vec = trim(strsplit(zonecodes,",")[[1]])
      vec = vec[!duplicated(vec) & vec!=""]
      div1DF[k,"ZoneCodes"] = paste(vec,collapse=",")
    }
  }
  write.table(div1DF, file="./data/ABS_data_by_DZN/DZN/div1_v3.csv", row.names = FALSE, sep="\t")

  print("==== all done.")
}

# given a zonecode, return all the job categories it can contain
findJobCatsByZoneCode <- function(zonecode){
  JobCats = c()
  #st = Sys.time()
  # for div 1
  vec = as.character(div1[,"ZoneCodes"])
  filter = rep(FALSE, length(vec))
  cnt = 1
  for(zcstr in vec){
    if(is.na(zcstr)){
      cnt = cnt + 1
      next
    }
    zcs = trim(unlist(strsplit(zcstr,",")))
    for(zc in zcs){
      matchrlt = regexpr(zc, zonecode, ignore.case=TRUE)
      if(attr(matchrlt,"match.length")==nchar(zonecode)){
        filter[cnt] = TRUE
        break
      }
    }
    cnt = cnt + 1
  }
  
  if(sum(filter) > 0){
    tmpJobCats = paste(div1[filter,"Div1Code"],"000",sep="")
    JobCats = c(JobCats, tmpJobCats)
  }
  
  # for div 2
  vec = as.character(div2[,"ZoneCodes"])
  filter = rep(FALSE, length(vec))
  cnt = 1
  for(zcstr in vec){
    if(is.na(zcstr)) {
      cnt = cnt + 1
      next
    }
    zcs = trim(unlist(strsplit(zcstr,",")))
    for(zc in zcs){
      matchrlt = regexpr(zc, zonecode, ignore.case=TRUE)
      if(attr(matchrlt,"match.length")==nchar(zonecode)){
        filter[cnt] = TRUE
        break
      }
    }
    cnt = cnt + 1
  }
  
  if(sum(filter) > 0){
    tmpJobCats = paste(div2[filter,"Div2Code"],"00",sep="")
    tmpJobCats = paste("X",tmpJobCats,sep="")
    JobCats = c(JobCats, tmpJobCats)
  }
  # for div 3
  vec = as.character(div3[,"ZoneCodes"])
  filter = rep(FALSE, length(vec))
  cnt = 1
  for(zcstr in vec){
    if(is.na(zcstr)) {
      cnt = cnt + 1
      next
    }
    zcs = trim(unlist(strsplit(zcstr,",")))
    for(zc in zcs){
      matchrlt = regexpr(zc, zonecode, ignore.case=TRUE)
      if(attr(matchrlt,"match.length")==nchar(zonecode)){
        filter[cnt] = TRUE
        break
      }
    }
    cnt = cnt + 1
  }
  
  if(sum(filter) > 0){
    tmpJobCats = paste(div3[filter,"Div3Code"],"0",sep="")
    tmpJobCats = paste("X",tmpJobCats,sep="")
    JobCats = c(JobCats, tmpJobCats)
  }
  # for div 4
  vec = as.character(div4[,"ZoneCodes"])
  filter = rep(FALSE, length(vec))
  cnt = 1
  for(zcstr in vec){
    if(is.na(zcstr)) {
      cnt = cnt + 1
      next
    }
    zcs = trim(unlist(strsplit(zcstr,",")))
    for(zc in zcs){
      matchrlt = regexpr(zc, zonecode, ignore.case=TRUE)
      if(attr(matchrlt,"match.length")==nchar(zonecode)){
        filter[cnt] = TRUE
        break
      }
    }
    cnt = cnt + 1
  }
  
  if(sum(filter) > 0){
    tmpJobCats = paste("X",div4[filter,"Div4Code"],sep="")
    JobCats = c(JobCats, tmpJobCats)
  }
  #et = Sys.time()
  #print(sprintf("==========> all done (in %.6f seconds) <==========", as.numeric(et-st, units="secs")))
  
  JobCats = JobCats[!duplicated(JobCats)]
  
  return(JobCats)
}

# given a job category, return the zonecodes that can contain it
findZoneCodesByJobCat <- function(jobcat){
  trimedJobCat = f_TrimJobCat(jobcat)
  divLevel = nchar(trimedJobCat)
  zonecodeStr = ""
  if(divLevel==1){
    vec = as.character(div1[,"Div1Code"])
    filter = grep(trimedJobCat, vec, ignore.case=TRUE, value=FALSE)
    zonecodeStr = trim(as.character(div1[filter,"ZoneCodes"]))
  } else if(divLevel==2){
    vec = as.character(div2[,"Div2Code"])
    filter = grep(trimedJobCat, vec, ignore.case=TRUE, value=FALSE)
    zonecodeStr = trim(as.character(div2[filter,"ZoneCodes"]))
  } else if(divLevel==3){
    vec = as.character(div3[,"Div3Code"])
    filter = grep(trimedJobCat, vec, ignore.case=TRUE, value=FALSE)
    zonecodeStr = trim(as.character(div3[filter,"ZoneCodes"]))
  }  else if(divLevel==4){
    vec = as.character(div4[,"Div4Code"])
    filter = grep(trimedJobCat, vec, ignore.case=TRUE, value=FALSE)
    zonecodeStr = trim(as.character(div4[filter,"ZoneCodes"]))
  }
  
  if(is.na(zonecodeStr) ) return(c())
  
  zonecodes = trim(unlist(strsplit(zonecodeStr,",")))
  return(zonecodes)
}

# trim job category string from "A000"->"A", "X0100"->"01", "X0110"->"011", "X0111"->"0111"
f_TrimJobCat <- function(jobcat){
  jobcatStr = jobcat
  jobcatStr = sub("X","",jobcatStr)
  jobcatStr = paste(jobcatStr, "#", sep="")
  if (regexpr("[A-W]000#",jobcatStr) > 0){
    return(sub("000#","",jobcatStr))
  }
  
  if (regexpr("[0-9][0-9]00#",jobcatStr) > 0){
    return(sub("00#","",jobcatStr))
  }
  
  if (regexpr("[0-9][0-9][0-9]0#",jobcatStr) > 0){
    return(sub("0#","",jobcatStr))
  }

  return(sub("#","",jobcatStr))
}

f_GenJobCatJSON <- function(){
  tmpStr = sprintf("{\"code\":\"%s000\", \"name\":\"%s\"}", div1[,"Div1Code"], div1[,"Div1Name"])
  write(sprintf("[%s]",paste(tmpStr, collapse=",")),"div1.txt")
  
  tmpStr = sprintf("{\"code\":\"X%s00\", \"name\":\"%s\", \"parent\":\"%s000\"}", div2[,"Div2Code"], div2[,"Div2Name"] , div2[,"Div1Code"])
  write(sprintf("[%s]",paste(tmpStr, collapse=",")),"div2.txt")
  
  tmpStr = sprintf("{\"code\":\"X%s0\", \"name\":\"%s\", \"parent\":\"X%s00\"}", div3[,"Div3Code"], div3[,"Div3Name"] , div3[,"Div2Code"])
  write(sprintf("[%s]",paste(tmpStr, collapse=",")),"div3.txt")
  
  tmpStr = sprintf("{\"code\":\"X%s\", \"name\":\"%s\", \"parent\":\"X%s0\"}", div4[,"Div4Code"], div4[,"Div4Name"] , div4[,"Div3Code"])
  write(sprintf("[%s]",paste(tmpStr, collapse=",")),"div4.txt")
}

f_GeoJSONWrite <- function(){
  
  outShp <- readShapePoly("./outputs/SplitPoly_X_Employment_fullcode")
  attr(outShp@proj4string,"projargs") = "+proj=longlat +ellps=GRS80 +no_defs"
  #CONST_EPSG4283_proj4string = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
  #outShp_pj = spTransform(outShp,CRS(CONST_EPSG4283_proj4string))
  outShp_pj = outShp
  outShp_pj@data = outShp_pj@data[,c("LGA_CODE","LGA","ZONE_CODE","X2310","X2412","X8500")]
  # generate KML output
  writeOGR(obj=outShp_pj, dsn="./outputs", layer="smalldata", driver = "ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)

  x<-readOGR(dsn="./outputs/SplitPoly_X_Employment_fullcode_selected.geojson", layer = 'OGRGeoJSON')
  CONST_projected_proj4string = "+proj=merc +datum=WGS84"
  x_pj = spTransform(x,CRS(CONST_projected_proj4string))
  writeOGR(obj=x_pj, dsn="./outputs/SplitPoly_X_Employment_fullcode_selected_pj.geojson", layer="tmp", driver = "GeoJSON", check_exists=TRUE, overwrite_layer=TRUE)
  
  # read geojson from url
  url = "./outputs/dump.geojson"
  layerFromUrl<-readOGR(dsn=url,layer = 'OGRGeoJSON')
  
  # read geojson from string
  jsonstr = readChar(url, file.info(url)$size)
  layerFromString<-readOGR(dsn=jsonstr,layer = 'OGRGeoJSON')
}

f_DataFrame2JSONString <- function(targetDataFrame=NULL)
{
  if(is.null(targetDataFrame)) return("[]")
  
  if(nrow(targetDataFrame)==0) return("[]")
  
  tmpFilePath = "./outputs/dump.geojson"
  tmpStr = ""
  # for spatial data frame, convert to geojson
  if(inherits(targetDataFrame, "Spatial")){
    if(file.exists(tmpFilePath)) file.remove(tmpFilePath)   
    writeOGR(obj=targetDataFrame, dsn=tmpFilePath, layer="tmp", driver = "GeoJSON",  check_exists=TRUE, overwrite_layer=TRUE)    
    tmpStr = readChar(tmpFilePath, file.info(tmpFilePath)$size)
    return(tmpStr)
  }
  
  # for normal data frame, convert to json
  vec=c()
  for(i in 1:nrow(targetDataFrame)){
    vec[i] = toJSON(targetDataFrame[i,], collapse="")
  }
  tmpStr = paste(vec,collapse=",")
  tmpStr = paste("[",tmpStr,"]",sep="")
  return(tmpStr)
  
}

f_GenBaseMap_2006 <- function(){
# algorithm start time
sTime = Sys.time()

# load shapefile in a much faster way 
dznempDF <- readShapePoly("./data/ABS_data_by_DZN/DZN/DZN_X_Employment_Orginal_TargetArea")
splitDF <- readShapePoly("./data/ABS_data_by_DZN/DZN/plan_zone_fullcode_split_checked_pvg2006_contain")
#dznempDF <- readOGR(dsn="./data/ABS_data_by_DZN/DZN/",layer="DZN_X_Employment_TargetArea",encoding="utf8")
#splitDF <- readOGR(dsn="./data/ABS_data_by_DZN/DZN/",layer="plan_zone_fullcode_split",encoding="utf8")

# create a set of shrunk polygons to ensure gContains can be performed accurately
shrinkedSplitPolyList = splitDF@polygons
for (j in 1:nrow(splitDF)){
  r = sqrt(splitDF[j,]@polygons[[1]]@area / pi)/10
  s = gBuffer(splitDF[j,], width = -r)@polygons[[1]]
  if (s@area==0){
    #print(sprintf("idx %i null after shrinking", j))
    shrinkedSplitPolyList[[j]]= splitDF[j,]@polygons[[1]]
  } else
  {
    shrinkedSplitPolyList[[j]]= s
  }
}
shrinkedSplitDF = splitDF
shrinkedSplitDF@polygons = shrinkedSplitPolyList

tmpSplitData = splitDF@data

# attach job categories to splitDF data
curDiv1Code = ""
curDiv2Code = ""
curDiv3Code = ""
curDiv4Code = ""
for (i in 1:nrow(div4)){
  if (div4[i, "Div1Code"] != curDiv1Code){
    curDiv1Code = div4[i, "Div1Code"]
    tmpSplitData[,paste(curDiv1Code,"000",sep="")] = 0
    #tmpSplitData = cbind(tmpSplitData, 0)
    #colnames(tmpSplitData)[length(colnames(tmpSplitData))] = paste(curDiv1Code,"000",sep="")
  }
  
  if (div4[i, "Div2Code"] != curDiv2Code){
    curDiv2Code = div4[i, "Div2Code"]
    newColName = paste("X",curDiv2Code,"00",sep="")
    if(!any(colnames(tmpSplitData)==newColName)){
      tmpSplitData[,newColName] = 0
    }
  }
  
  if (div4[i, "Div3Code"] != curDiv3Code){
    curDiv3Code = div4[i, "Div3Code"]
    newColName = paste("X",curDiv3Code,"0",sep="")
    if(!any(colnames(tmpSplitData)==newColName)){
      tmpSplitData[,newColName] = 0
    }
  }
  
  if (div4[i, "Div4Code"] != curDiv4Code){
    curDiv4Code = div4[i, "Div4Code"]
    newColName = paste("X",curDiv4Code,sep="")
    if(!any(colnames(tmpSplitData)==newColName)){
      tmpSplitData[,newColName] = 0
    }
  }
}

# convert into matrix for speed
tmpSplitDataMat = as.matrix(tmpSplitData)

for (i in 1:nrow(dznempDF)){ 
  print(sprintf("=== processing: %i of %i",i, nrow(dznempDF)))
  dznArea = dznempDF[i,]
  containedIdx = c()
  
  x = shrinkedSplitDF %over% dznArea
  names(x) = c(1:length(x))
  subx = x[!is.na(x)]
  containedIdx = as.integer(names(subx))
  
  if (length(containedIdx) == 0) {
    print("------ no contained polygons found, skipped --")
    next
  }
  
  # get non-duplicated zonecodes within that dnzarea
  zonecodes = levels(as.factor(tmpSplitDataMat[containedIdx,"ZONE_CODE"]))
  
  # find jobcats containing these zonecodes
  JobCats = c()
  for (zc in zonecodes){
    
    tmpJobCats = findJobCatsByZoneCode(zc)
    if(is.null(tmpJobCats)) next
    
    JobCats = c(JobCats, tmpJobCats)
  }
  
  if (is.null(JobCats)) next
  
  JobCats = JobCats[!duplicated(JobCats)]
  
  # only interest in these jobcats which contain jobs
  filter = as.vector(dznempDF@data[i,JobCats]>0)
  JobCats = JobCats[filter]
  
  if (length(JobCats) == 0) next
  
  # assign jobs in each jc to split polygons based on zonecode matching
  for (jc in JobCats){
    # find zodecode(s) that jc contains
    zcs = findZoneCodesByJobCat(jc)
  
    if (length(zcs) == 0) next
    
    subContainedIdx = c()
    subContainedArea = c()
    
    for (cidx in containedIdx){      
      for (zc in zcs){
        matchrlt = regexpr(zc, tmpSplitDataMat[cidx,"ZONE_CODE"], ignore.case=TRUE)
        if (attr(matchrlt,"match.length")==nchar(tmpSplitDataMat[cidx,"ZONE_CODE"])){
          subContainedIdx = c(subContainedIdx, cidx)
          subContainedArea = c(subContainedArea, splitDF@polygons[[cidx]]@area)
          break
        }
      }
    }
    
    subSumArea = sum(subContainedArea)
    subContainedJobNumsRatio = subContainedArea / subSumArea
    for (j in 1:length(subContainedIdx)){
        tmpSplitDataMat[subContainedIdx[j], jc] = round(as.numeric(dznempDF@data[i,jc]) * subContainedJobNumsRatio[j], digits = 4)
    }
  }

}



tmpmx = as.data.frame(tmpSplitDataMat)
# collapse jobnumbers
tmpmx[which(is.na(tmpmx))] = 0

cols = colnames(tmpmx)
numericCols = cols[!(c(1:length(cols)) %in% c(4,6,8))] # c(3,5,7,8,10)
for(i in 1:length(numericCols)){
  print(sprintf("coverting numeric columns %i of %i", i, length(numericCols)))
  tmpmx[,numericCols[i]] = as.numeric(as.vector(tmpmx[,numericCols[i]]))
}


# sum up on the 3 level
colnames = paste("X",div3[,"Div3Code"],"0", sep="")
divcodes = div3[,"Div3Code"]
for(i in 1:length(colnames)){
  colname = colnames[i]
  divcode = divcodes[i]
  filter = div4[,"Div3Code"] == divcode
  subcolnames = paste("X",div4[filter,"Div4Code"], sep="")
  if(subcolnames != colname){
    tmpmx[,colname] = tmpmx[,colname] + rowSums(tmpmx[,subcolnames, drop=FALSE])
  }
}

# sum up on the 2 level
colnames = paste("X",div2[,"Div2Code"],"00", sep="")
divcodes = div2[,"Div2Code"]
for(i in 1:length(colnames)){
  colname = colnames[i]
  divcode = divcodes[i]
  filter = div3[,"Div2Code"] == divcode
  subcolnames = paste("X",div3[filter,"Div3Code"],"0", sep="")
  if(subcolnames != colname){
    tmpmx[,colname] = tmpmx[,colname] + rowSums(tmpmx[,subcolnames, drop=FALSE])
  }
}

# sum up on the 1 level
colnames = paste(div1[,"Div1Code"],"000", sep="")
divcodes = div1[,"Div1Code"]
for(i in 1:length(colnames)){
  colname = colnames[i]
  divcode = divcodes[i]
  filter = div2[,"Div1Code"] == divcode
  subcolnames = paste("X",div2[filter,"Div2Code"],"00", sep="")
  tmpmx[,colname] = tmpmx[,colname] + rowSums(tmpmx[,subcolnames, drop=FALSE])
}


# save the result
splitDF@data = tmpmx

writeOGR(obj=splitDF, dsn="./outputs", layer="MidPoly_X_Employment_2006", driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)

# algorithm end time
eTime = Sys.time()

print(sprintf("==========> all done (in %.2f seconds) <==========", as.numeric(eTime-sTime, units="secs")))
}

f_GenBaseMap_2011 <- function(){
  # algorithm start time
  sTime = Sys.time()
  
  # load shapefile in a much faster way 
  dznempDF <- readShapePoly("./data/ABS_data_by_DZN_2011/DZN_X_Employment_2011_Orgina_TargetArea")
  splitDF <- readShapePoly("./data/ABS_data_by_DZN_2011/plan_zone_fullcode_split_2011_checked_pvg2011_contain")
  #dznempDF <- readOGR(dsn="./data/ABS_data_by_DZN/DZN/",layer="DZN_X_Employment_TargetArea",encoding="utf8")
  #splitDF <- readOGR(dsn="./data/ABS_data_by_DZN/DZN/",layer="plan_zone_fullcode_split",encoding="utf8")
  
  # create a set of shrunk polygons to ensure gContains can be performed accurately
  shrinkedSplitPolyList = splitDF@polygons
  for (j in 1:nrow(splitDF)){
    r = sqrt(splitDF[j,]@polygons[[1]]@area / pi)/10
    s = gBuffer(splitDF[j,], width = -r)@polygons[[1]]
    if (s@area==0){
      #print(sprintf("idx %i null after shrinking", j))
      shrinkedSplitPolyList[[j]]= splitDF[j,]@polygons[[1]]
    } else
    {
      shrinkedSplitPolyList[[j]]= s
    }
  }
  shrinkedSplitDF = splitDF
  shrinkedSplitDF@polygons = shrinkedSplitPolyList
  
  tmpSplitData = splitDF@data
  
  # attach job categories to splitDF data
  curDiv1Code = ""
  curDiv2Code = ""
  curDiv3Code = ""
  curDiv4Code = ""
  for (i in 1:nrow(div4)){
    if (div4[i, "Div1Code"] != curDiv1Code){
      curDiv1Code = div4[i, "Div1Code"]
      tmpSplitData[,paste(curDiv1Code,"000",sep="")] = 0
      #tmpSplitData = cbind(tmpSplitData, 0)
      #colnames(tmpSplitData)[length(colnames(tmpSplitData))] = paste(curDiv1Code,"000",sep="")
    }
    
    if (div4[i, "Div2Code"] != curDiv2Code){
      curDiv2Code = div4[i, "Div2Code"]
      newColName = paste("X",curDiv2Code,"00",sep="")
      if(!any(colnames(tmpSplitData)==newColName)){
        tmpSplitData[,newColName] = 0
      }
    }
    
    if (div4[i, "Div3Code"] != curDiv3Code){
      curDiv3Code = div4[i, "Div3Code"]
      newColName = paste("X",curDiv3Code,"0",sep="")
      if(!any(colnames(tmpSplitData)==newColName)){
        tmpSplitData[,newColName] = 0
      }
    }
    
    if (div4[i, "Div4Code"] != curDiv4Code){
      curDiv4Code = div4[i, "Div4Code"]
      newColName = paste("X",curDiv4Code,sep="")
      if(!any(colnames(tmpSplitData)==newColName)){
        tmpSplitData[,newColName] = 0
      }
    }
  }
  
  # convert into matrix for speed
  tmpSplitDataMat = as.matrix(tmpSplitData)
  
  for (i in 1:nrow(dznempDF)){ 
    print(sprintf("=== processing: %i of %i",i, nrow(dznempDF)))
    dznArea = dznempDF[i,]
    containedIdx = c()
    
    x = shrinkedSplitDF %over% dznArea
    names(x) = c(1:length(x))
    subx = x[!is.na(x)]
    containedIdx = as.integer(names(subx))
    
    if (length(containedIdx) == 0) {
      print("------ no contained polygons found, skipped --")
      next
    }
    
    # get non-duplicated zonecodes within that dnzarea
    zonecodes = levels(as.factor(tmpSplitDataMat[containedIdx,"ZONE_CODE"]))
    
    # find jobcats containing these zonecodes
    JobCats = c()
    for (zc in zonecodes){
      
      tmpJobCats = findJobCatsByZoneCode(zc)
      if(is.null(tmpJobCats)) next
      
      JobCats = c(JobCats, tmpJobCats)
    }
    
    if (is.null(JobCats)) next
    
    JobCats = JobCats[!duplicated(JobCats)]
    
    # only interest in these jobcats which contain jobs
    filter = as.vector(dznempDF@data[i,JobCats]>0)
    JobCats = JobCats[filter]
    
    if (length(JobCats) == 0) next
    
    # assign jobs in each jc to split polygons based on zonecode matching
    for (jc in JobCats){
      # find zodecode(s) that jc contains
      zcs = findZoneCodesByJobCat(jc)
      
      if (length(zcs) == 0) next
      
      subContainedIdx = c()
      subContainedArea = c()
      
      for (cidx in containedIdx){      
        for (zc in zcs){
          matchrlt = regexpr(zc, tmpSplitDataMat[cidx,"ZONE_CODE"], ignore.case=TRUE)
          if (attr(matchrlt,"match.length")==nchar(tmpSplitDataMat[cidx,"ZONE_CODE"])){
            subContainedIdx = c(subContainedIdx, cidx)
            subContainedArea = c(subContainedArea, splitDF@polygons[[cidx]]@area)
            break
          }
        }
      }
      
      subSumArea = sum(subContainedArea)
      subContainedJobNumsRatio = subContainedArea / subSumArea
      for (j in 1:length(subContainedIdx)){
        tmpSplitDataMat[subContainedIdx[j], jc] = round(as.numeric(dznempDF@data[i,jc]) * subContainedJobNumsRatio[j], digits = 4)
      }
    }
    
  }
  
  
  
  tmpmx = as.data.frame(tmpSplitDataMat)
  # collapse jobnumbers
  tmpmx[which(is.na(tmpmx))] = 0
  
  cols = colnames(tmpmx)
  numericCols = cols[!(c(1:length(cols)) %in% c(4,6,8))] # c(3,5,7,8,10)
  for(i in 1:length(numericCols)){
    print(sprintf("coverting numeric columns %i of %i", i, length(numericCols)))
    tmpmx[,numericCols[i]] = as.numeric(as.vector(tmpmx[,numericCols[i]]))
  }
  
  
  # sum up on the 3 level
  colnames = paste("X",div3[,"Div3Code"],"0", sep="")
  divcodes = div3[,"Div3Code"]
  for(i in 1:length(colnames)){
    colname = colnames[i]
    divcode = divcodes[i]
    filter = div4[,"Div3Code"] == divcode
    subcolnames = paste("X",div4[filter,"Div4Code"], sep="")
    if(subcolnames != colname){
      tmpmx[,colname] = tmpmx[,colname] + rowSums(tmpmx[,subcolnames, drop=FALSE])
    }
  }
  
  # sum up on the 2 level
  colnames = paste("X",div2[,"Div2Code"],"00", sep="")
  divcodes = div2[,"Div2Code"]
  for(i in 1:length(colnames)){
    colname = colnames[i]
    divcode = divcodes[i]
    filter = div3[,"Div2Code"] == divcode
    subcolnames = paste("X",div3[filter,"Div3Code"],"0", sep="")
    if(subcolnames != colname){
      tmpmx[,colname] = tmpmx[,colname] + rowSums(tmpmx[,subcolnames, drop=FALSE])
    }
  }
  
  # sum up on the 1 level
  colnames = paste(div1[,"Div1Code"],"000", sep="")
  divcodes = div1[,"Div1Code"]
  for(i in 1:length(colnames)){
    colname = colnames[i]
    divcode = divcodes[i]
    filter = div2[,"Div1Code"] == divcode
    subcolnames = paste("X",div2[filter,"Div2Code"],"00", sep="")
    tmpmx[,colname] = tmpmx[,colname] + rowSums(tmpmx[,subcolnames, drop=FALSE])
  }
  
  
  # save the result
  splitDF@data = tmpmx
  
  writeOGR(obj=splitDF, dsn="./outputs", layer="MidPoly_X_Employment_2011", driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
  
  # algorithm end time
  eTime = Sys.time()
  
  print(sprintf("==========> all done (in %.2f seconds) <==========", as.numeric(eTime-sTime, units="secs")))
}


