### BEGIN EXAMPLE aggregarting polygons using dissolve --------------------------------------------------------------------
# from: http://www.nceas.ucsb.edu/scicomp/usecases/PolygonDissolveOperationsR

# required packages
library(maptools)   # for geospatial services; also loads foreign and sp
library(gpclib)     # General Polygon Clipping library 
library(rgdal)      # for map projection work; also loads sp
library(rgeos)
library(PBSmapping) # for GIS_like geospatial object manipulation / anslysis including poly
gpclibPermit()
require(gpclib)

f_norm <- function(vec, na.rep.auto = FALSE, na.rep = NA){
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # normalize a numeric vector, if na.rep.auto = TRUE, the normalized value for NA elements will be calculated using: (mean(vec) - min(vec))/(max(vec) - min(vec)) and the na.rep will be ignored
  # otherwise, na.rep (clap to 0-1) will be used as the normalized value for NA elements.
  # usage:
  # (1) f_norm(x):
  # (2) f_norm(x, na.rep=0.667)
  # (3) f_norm(x, na.rep.auto=TRUE)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # find the max and min, NA elements are removed (otherwise, NA will be returned if x contains NA elements)
  
  # do some inputs type tests before starting
  if (mode(vec) != "numeric"){
    stop("a numeric vector is expected for vec")
  }
  
  if (mode(na.rep.auto) != "logical" || length(na.rep.auto) != 1){
    na.rep.auto = FALSE
    warning("TRUE/FALSE is expected for na.rep.auto. default value FALSE is applied")
  }
  
  if (!is.na(na.rep) && (mode(na.rep) != "numeric" || length(na.rep) != 1)){
    na.rep = NA
    warning("a numeric number or NA is expected for na.rep. default value NA is applied")
  }
  
  v_max = max(vec, na.rm = TRUE)
  v_min = min(vec, na.rm = TRUE)
  
  # if na.rep is assigned, make sure it is between 0-1
  if(!is.na(na.rep) && (na.rep < 0 || na.rep >1)){
    na.rep = 0.5
  }
  
  # if vec is empty or all elements in vec are NA, use na.rep as the value of the normalized output
  if (v_max == -Inf || v_min == -Inf){
    return(rep(na.rep,length(vec)))
  }
  
  # if max == min, set normalized value to 1 for all vector elements
  if (v_max == v_min){
    rtn = vec - v_min + 1
  } else
  {
    # do the normalization
    rtn = (vec - v_min)/(v_max - v_min)
  }
  
  # if na.rep.auto, try to use normalized the mean and use the value as na.rep for output
  if (na.rep.auto){
    if (v_max == v_min){
      na.rep = 1
    } else{
      na.rep = (mean(vec[!is.na(vec)]) - v_min) / (v_max - v_min)
    }
  }
  
  # replace any NA with na.rep
  rtn[is.na(rtn)] = na.rep
  
  return(rtn)
}


f_wards <- function(adata, pdata, ianmwh, snswh=c(0.5,0.5)){
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # do ward's clustering on spatial and non-spatial attributes.
  # input:
  # (1) adata: attribute data (data.frame)
  # (2) pdata: polygons (list)
  # (3) ianmwh: interested attribute data column names and wreights (data.frame) e.g.d = data.frame(cbind(ATTR_NAME=c("Morans_I","attr1","attr2"), ATTR_WEIGHT=c(0.4,0.3,0.3)))
  # (4) snswh: spatial and non-spatial weights (vector)
  # output:
  # new shp file with a new data column marks its ward's cluster number
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # step1. inputs validation
  
  # step2. backup data, adata and pdata will be changed during clustering
  adata_bak = adata
  pdata_bak = pdata
  
  IN_NS_ATTR_NAMES = as.character(ianmwh[[1]])
  IN_NS_ATTR_WEIGHTS = as.numeric(as.character(ianmwh[[2]]))
  # step3. do the ward's clustering
  while(length(adata[[1]]) > 1) #continue if more than one rows exist
  {
    # reset these variables each loop
    min_pdf_dist = 999999999
    min_idx_i = -1
    min_idx_j = -1
    
    CUR_POLYGON_NUM = length(adata[[1]])
    print(sprintf("=============================== merging, %d plogyons remain",CUR_POLYGON_NUM))                         
    # normalize attribute data
    # find the interested non-spatial attrs
    norm_ns_attrs_dataframe = adata[IN_NS_ATTR_NAMES]
    # then do the normalization, NA will be retained in the result
    for(ns_attr_no in 1:length(norm_ns_attrs_dataframe)){
      norm_ns_attrs_dataframe[[ns_attr_no]] = f_norm(norm_ns_attrs_dataframe[[ns_attr_no]])
    }
    
    # init a dataframe to record non-normalized s_dist and ns_dist, as well as the polygon idx pair: idx_i, idx_j  
    len_dist_dataframe = CUR_POLYGON_NUM * (CUR_POLYGON_NUM - 1) / 2
    dist_dataframe = data.frame(cbind(s_dist=rep(0.5,len_dist_dataframe), ns_dist=rep(0.5,len_dist_dataframe), idx_i=rep(0,len_dist_dataframe),  idx_j=rep(0,len_dist_dataframe)))
    dist_dataframe_rowcounter = 1
    for (idx_i in 1:(CUR_POLYGON_NUM-1)){
      
      pi = SpatialPolygons(list(pdata[[idx_i]]))
      len_pi = gLength(pi)
      
      for (idx_j in (idx_i+1):CUR_POLYGON_NUM){
        
        # calc spatial distance
        pj = SpatialPolygons(list(pdata[[idx_j]]))
        len_pj = gLength(pj)
        hdist =  gDistance(pi,pj,hausdorff=TRUE)
        
        pij = SpatialPolygons(list(pdata[[idx_i]],pdata[[idx_j]]))
        union_pij = unionSpatialPolygons(pij ,c(1,1))
        len_union_pij = gLength(union_pij)
        len_shared_boundary_pij = len_pi + len_pj - len_union_pij
        
        # ref: Joshi' paper p. 22
        sp_dist = hdist * (1 - 2 * len_shared_boundary_pij / (len_pi + len_pj))
        
        print(sprintf("(p%d,p%d) spatial distance -- hdist_adjust: %.2f,    shared: %.2f,    hdist: %.2f",idx_i, idx_j, sp_dist, len_shared_boundary_pij, hdist))
        
        # calc non spatial distance (Euclidian distance)
        # pre-conditions: (1) attribute is numeric (2) attribute is normalized (i.e., 0-1)
        # issue: how to handle the distance between NA attribute value(s)? A simple way is to set the distance 0 so it makes no contribution (not exactly, it actually shortens the "real" distance) to the final dPDF.
        
        nsp_dist = 0
        for(ns_attr_no in 1:length(norm_ns_attrs_dataframe)){
          nsp_dist = nsp_dist + IN_NS_ATTR_WEIGHTS[ns_attr_no]*(norm_ns_attrs_dataframe[idx_i,][ns_attr_no][[1]] - norm_ns_attrs_dataframe[idx_j,][ns_attr_no][[1]])^2
        }
        nsp_dist = sqrt(nsp_dist)
        
        # set row values for dist_dataframe
        dist_dataframe[dist_dataframe_rowcounter,]["s_dist"][[1]] = sp_dist
        dist_dataframe[dist_dataframe_rowcounter,]["ns_dist"][[1]] = nsp_dist
        dist_dataframe[dist_dataframe_rowcounter,]["idx_i"][[1]] = idx_i
        dist_dataframe[dist_dataframe_rowcounter,]["idx_j"][[1]] = idx_j
        
        dist_dataframe_rowcounter = dist_dataframe_rowcounter + 1
        
        print(sprintf("(p%d,p%d) non-spatial distance --  %.6f",idx_i, idx_j, nsp_dist))
        
      }
    }
    
    # then let's normalize the s_dist (ns_dist) before computing the dPDF
    #dist_dataframe[["s_dist"]] = f_norm(dist_dataframe[["s_dist"]])
    # the meaning of normalized distance value on 'normalized attribtes' is vague.  need further disccussion if I cannot figure out
    #dist_dataframe[["ns_dist"]] = f_norm(dist_dataframe[["ns_dist"]])
    
    # now we can compute the dPDF and find the smallest one
    for (dist_idx in 1:len_dist_dataframe){
      pdf_dist = snswh[1]*dist_dataframe[dist_idx,]["s_dist"][[1]] + snswh[2]*dist_dataframe[dist_idx,]["ns_dist"][[1]]
      if (pdf_dist < min_pdf_dist){
        min_pdf_dist = pdf_dist
        min_idx_i = dist_dataframe[dist_idx,]["idx_i"][[1]]
        min_idx_j = dist_dataframe[dist_idx,]["idx_j"][[1]]
      }
    }
    
    print(sprintf("min pdf distance found between (p%d,p%d) :  %.6f",min_idx_i, min_idx_j,min_pdf_dist))
    
    if (min_pdf_dist > CONST_DISTANCE_THRESHOLD){
      # cannot aggreate any more, should exit algorithm
      print("********** merge condition not meet, exit.")
      break
    } else {
      # merge two polygons: p(min_idx_i) and p(min_idx_j)
      
      # step1. handle the geometry data
      tmp_pij = SpatialPolygons(list(pdata[[min_idx_i]],pdata[[min_idx_j]]))
      uni_pij = unionSpatialPolygons(tmp_pij ,c(1,1))
      uni_pij@polygons[[1]]@ID = as.character(min_idx_i)
      
      # store the newly merged polygon to DissolvePolyList[[min_idx_i]]
      pdata[min_idx_i] = uni_pij@polygons
      # then remove DissolvePolyList[[min_idx_j]] from the list
      pdata[min_idx_j] = NULL
      
      # step2. handle the attribute data
      # merge attribute data into min_idx_i
      ns_attr_no = 1
      for(ns_attr_no in 1:length(norm_ns_attrs_dataframe)){
        adata[min_idx_i,][IN_NS_ATTR_NAMES[ns_attr_no]][[1]] = adata[min_idx_i,][IN_NS_ATTR_NAMES[ns_attr_no]][[1]] +  adata[min_idx_j,][IN_NS_ATTR_NAMES[ns_attr_no]][[1]]
      }
      
      # find the wardclust value of row min_idx_j in dataTargetYear before remove it
      wardclust_j_val = adata[min_idx_j,]["wardclut"][[1]]
      # find the wardclust value that will be used to replace wardclust_j_val in the dataTargetYear_bak
      wardclust_i_val = adata[min_idx_i,]["wardclut"][[1]]
      
      # remove data row on min_idx_j
      adata = adata[c(-min_idx_j),]
      
      # update wardclust column in the dataTargetYear_bak
      # find all rows in dataTargetYear_bak whose wardclust value == wardclust_j_val
      tmpfilter = adata_bak[,]["wardclut"][[1]] == wardclust_j_val
      adata_bak[tmpfilter,]["wardclut"][[1]] = wardclust_i_val
      
      # visualize results of each clustering step
      # should be commented when released
      # plot(pdata) 
    }
  }
}

#set up the working directory, ATTENTION: the data folder used below is "District497", no whitespaces
setwd("/Users/Shared/Documents/AURIN/R")
data = read.table("./data/Modified Data Tables/Morans_I/Morans_I_type1_urbanonly.txt", header=TRUE, na.strings="NA", blank.lines.skip = TRUE, fill = TRUE)

# constant variables definition
CONST_MISSING_DATA_DISTANCE = 0.5
CONST_LARGE_CLUSTER_IDX = 9999



# using a year-filter to select the data we want to analyze
targetYear = 2001
targetYearFilterVec = data[[1]] == targetYear
filteredDataTargetYear = data[targetYearFilterVec,]
# create cluster idx columns to hold cluster result
clustIdxCol = c(10,20,30,40,50,60,70,80)
clustIdxColNames = paste("C",clustIdxCol,sep="")
for (i in 1:length(clustIdxColNames)){
  filteredDataTargetYear[,clustIdxColNames[i]]=0
}

# get all valid data row by removing those containing NA value
dataTargetYear = filteredDataTargetYear[complete.cases(filteredDataTargetYear),]
# backup the NA data row if exists
NA_dataTargetYear = filteredDataTargetYear[!complete.cases(filteredDataTargetYear),]

# print the new data structure
#attr(dataTargetYear,"names")

# calc the distance on column "Morans_I" using "euclidean" distance
d = dist(x=dataTargetYear[[3]], method="euclidean")

# cluster d using ward linkage criterion. ATTENTION, d^2 is used for 'ward' method; for 'single' or 'complete' method, using d (not squared).
newtree <- hclust(d^2, method="ward")

# cut tree into k clusters, and store the clusterIdx
for (i in 1:length(clustIdxColNames)){
  dataTargetYear[,clustIdxColNames[i]] = cutree(newtree, k=clustIdxCol[i])
}


# now let's take care of the NA-rows
NA_ROW_NUN = dim(NA_dataTargetYear)[1]

if (NA_ROW_NUN>0){
  # mark each NA-row as an individual cluseter
  for (i in 1:length(clustIdxColNames)){
    NA_dataTargetYear[,clustIdxColNames[i]] = c((CONST_LARGE_CLUSTER_IDX-NA_ROW_NUN+1):CONST_LARGE_CLUSTER_IDX)    
  }
  
  # bind the NA_dataTargetYear back into dataTargetYear
  dataTargetYear = rbind(dataTargetYear,NA_dataTargetYear)
}

# add an additional column to record ward's cluster results
# by default, each polygon is a cluster
dataTargetYear[,"wardclut"]= c(1:length(dataTargetYear[[1]]))
# add two test attributes rows for 
# attr1 is random number (fractional) between 0.0 to 1.0
dataTargetYear[,"attr1"]= runif(length(dataTargetYear[[1]]), 0.0,1.0)
# attr2 is random number (integer) between 1 to 50, repeats are allowed
dataTargetYear[,"attr2"]= sample(1:50, length(dataTargetYear[[1]]), replace =T)

# make a backup here
dataTargetYear_bak = dataTargetYear

# Using the id_497 polygons for Indonesia, generate 
print("Reading and transforming id_497 shapefile into id_292...")
districts_idm <- readOGR(dsn="./data/shapefiles/District497",layer="crosswalk_all497_feb11",input_field_name_encoding="utf8")
original_proj4string = attr(districts_idm@proj4string,"projargs")
CONST_projected_proj4string = "+proj=merc +datum=WGS84"
# check if the original data is projected
# Transform the polygons (which were read in as unprojected geographic coordinates) to an Albers Equal Area projection
# this is essential since we need calculate distance/area for polygons
if (is.projected(districts_idm)){
  districts_idm_pj = districts_idm
} else
{
  # CRS(coordinate reference system http://itia.ntua.gr/antonis/technical/coordinate-systems)
  #
  # In our case, the ONLY purpose for CRS conversion, from GeodeticCRS to a ProjectedCRS, is to calculate geometry features (e.g., polyline areas, perimeter, shared boundary length) for difference matrix computation. The projected map CANNOT be used for visualization in the AURIN platform.
  # So maybe it's unnecessary (or difficult or impossible) to find a "proper" (or specific) proj4 strings (i.e, a set of projection parameters) for each input dataset
  # To keep simple, we use the most basic form ("+proj=merc +datum=WGS84") of proj4 string to do the conversion:
  districts_idm_pj = spTransform(districts_idm,CRS(CONST_projected_proj4string))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Here is the way to make a formal CRS conversion:
  # 1. find a "proper" proj4 strings for conversion, go to http://www.epsg-registry.org and search EPGS code by Area name, e.g. Australia, Indonesia
  # 2. pick a "proper" EPGS code from the search result list, e.g., 3001
  # 3. in the R console, print(CRSargs(CRS("+init=epsg:3001")))  (http://127.0.0.1:12201/help/library/sp/html/CRS-class.html) to get the proj4 string we need. e.g., "+proj=merc +lon_0=110 +k=0.997 +x_0=3900000 +y_0=900000 +ellps=bessel +towgs84=-377,681,-50,0,0,0,0 +units=m +no_defs"
  # 4. do the conversion using 'spTransform' method (in 'rgdal' package) with proj4 string:
  #   idm_pj = spTransform(idm,CRS("+proj=merc +lon_0=110 +k=0.997 +x_0=3900000 +y_0=900000 +ellps=bessel +towgs84=-377,681,-50,0,0,0,0 +units=m +no_defs"))
  # # # # # # # # # # # # # # # # # # # # # # # # # # #
}


# extract "parent_292" to a vector
col_292_vec = districts_idm_pj@data[, 34]

#get the groups for the vector, and convert to integer
col_292_vec_group = as.integer(levels(factor(col_292_vec)))

#create a list to store either dissolved or original polygons object
DissolvePolyList = list()

# number of raw data row length(col_292_vec_group)
TEST_DATA_ROW_NUM = 5
for (iCtr in 1:TEST_DATA_ROW_NUM){
  print(sprintf("%d==== is processing",iCtr))
  
  # for each loop, create a logical vector for dissolve use
  # if items in the logical vector are marked as "TRUE", their cooresponding polygon will be dissolved 
  select = districts_idm_pj@data[, 34] == col_292_vec_group[iCtr]
  
  # use the 'select' logical vector as a filter on the  districts_idm data rows
  groupedPolys = districts_idm_pj[select,]
  
  # if more than one rows are selected, dissolve function will get called
  if (length(groupedPolys) > 1){
    lps <- getSpPPolygonsLabptSlots(groupedPolys)
    IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
    Dissolve <- unionSpatialPolygons(groupedPolys ,IDOneBin)
    
    Dissolve@polygons[[1]]@ID = as.character(iCtr)
    DissolvePolyList[iCtr] = Dissolve@polygons
    
  } else #otherwise, append original one to the DissolvePolyVec
  {
    groupedPolys@polygons[[1]]@ID = as.character(iCtr)
    DissolvePolyList[iCtr] = groupedPolys@polygons
  }   
}

# make a backup here
DissolvePolyList_bak = DissolvePolyList

# build the geometry features (i.e., polygons)
Sr =SpatialPolygons(DissolvePolyList)

# setup a proper proj4string for outputs
if (is.projected(districts_idm)){
  Sr@proj4string = CRS(original_proj4string)
} else
{
  Sr@proj4string = CRS(CONST_projected_proj4string)
}

# create a new feature class holding both geometry features and attribute data
newDataFrame = SpatialPolygonsDataFrame(Sr,data=dataTargetYear[1:TEST_DATA_ROW_NUM,], match.ID = FALSE)


# draw the new map, the visualization part, convert back to original CRS if necessary
if (is.projected(districts_idm)){
  newDataFrame_original_proj = newDataFrame
} else
{
  newDataFrame_original_proj = spTransform(newDataFrame,CRS(original_proj4string))
}

newDataFrame_PS = SpatialPolygons2PolySet(newDataFrame_original_proj)
for (i in 1:length(clustIdxColNames))
{
  print(paste("producing map on with cluster number:",as.character(clustIdxCol[i])))
  plotPolys(newDataFrame_PS, proj = TRUE,col=dataTargetYear[[clustIdxColNames[i]]][1:TEST_DATA_ROW_NUM],xlab=paste("longitude",as.character(clustIdxCol[i])),ylab="latitude")  
}


# maybe we need change back to the original projection or 
# output data
writeOGR(obj=newDataFrame, dsn="./outputs", layer="newDataFrame", driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)

print("========= done.")

#################################
##  new program start

CONST_DISTANCE_THRESHOLD = 999.0 

# for test
dataTargetYear = dataTargetYear[1:TEST_DATA_ROW_NUM,]
DissolvePolyList = DissolvePolyList[1:TEST_DATA_ROW_NUM]
nmwt=data.frame(cbind(ATTR_NAME=c("Morans_I","attr1","attr2"), ATTR_WEIGHT=c(0.4,0.3,0.3)))

f_wards(dataTargetYear, DissolvePolyList,ianmwh = nmwt, snswh=c(0.005,0.995))
print("=========done")


