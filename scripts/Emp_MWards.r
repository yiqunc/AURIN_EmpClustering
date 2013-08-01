# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: the algorithm (f_wards) is developed to bulid clusters on conditions and visualise the results on map. 
# Version: 2.2
# Last Updated: 12-June-2013
# Supervised by: Jennifer Day    jennieday@gmail.com
# Coded by: Yiqun(Benny) Chen    chen_yiqun@hotmail.com
#
# change logs: v2.2
# (1) case insentivte of column names
# (2) change VC mode processing logic
# change logs: v2.1
# (1) add gDataYear to switch between datasets in different years, in batch mode, the output folder structure is : output/YYYY/dirname/gdist/spatialweight
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# required packages
library(maptools)   # for geospatial services; also loads foreign and sp
library(gpclib)     # General Polygon Clipping library 
library(rgdal)      # for map projection work; also loads sp
library(rgeos)
library(PBSmapping) # for GIS_like geospatial object manipulation / anslysis including poly
library(foreign) # for stata format output
library(RColorBrewer) # for color rendering on map 
gpclibPermit()
require(gpclib)
CONST_projected_proj4string = "+proj=merc +datum=WGS84"
# the output projection string is EPSG4283, which can be obtained this way:
# EPSG = make_EPSG()
# filter = EPSG[,"code"]=="4283"
# filter[which(is.na(filter))] = FALSE
# projstring = EPSG[filter,"prj4"]
CONST_EPSG4283_proj4string = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"

CONST_na_nsp_distance = 0.5
GLOBAL_polygon_id_prefix = "ANDS_M_"
GLOBAL_polygon_id_counter  = 1

#set up the working directory
setwd("/Users/Shared/Documents/AURIN/R")

# create a map color vector for rendering, 55 valid colors
gColVec=c(brewer.pal(11, "Spectral"), 
         brewer.pal(12, "Set3"),
         brewer.pal(8, "Set2"),
         brewer.pal(9, "Set1"),
         brewer.pal(8, "Pastel2"),
         brewer.pal(9, "Pastel1")
)

gColVec = gColVec[!duplicated(gColVec)] #55 different colors
gColVec = rep(gColVec,10) # repeat 10 times to generate a big enough number of colors
gReservedColVec = c("#B22222","#228B22","#FF8C00","#0000CD","#FF1493","#008000","#00FFFF","#0000FF","#808080","#000080","#FF00FF")
gColVec = c(gReservedColVec,gColVec)
gOriginal_proj4string = ""
gX = NULL
gX_pj = NULL

gDataYear = "2006"
gShowDebugInfo = TRUE
debugPrint <- function(str){
  if(gShowDebugInfo){
    print(str)
  }
}

f_getCentroidDistance <- function(latlon_s, latlon_e){
  return(sqrt((latlon_s[1]-latlon_e[1])^2 + (latlon_s[2]-latlon_e[2])^2))
} 

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
  } else {
    # do the normalization
    rtn = (vec - v_min)/(v_max - v_min)
  }
  
  # if na.rep.auto, try to use normalized the mean and use the value as na.rep for output
  if (na.rep.auto){
    if (v_max == v_min){
      na.rep = 1
    } else {
      na.rep = (mean(vec[!is.na(vec)]) - v_min) / (v_max - v_min)
    }
  }
  
  # replace any NA with na.rep
  rtn[is.na(rtn)] = na.rep
  
  return(rtn)
}

# wards: applied geodistance threshold, with distance calculation process optimised
f_wards <- function(adata, pdata, ianmwh, snswh=c(0.5,0.5), dthresh, proj4string=CONST_projected_proj4string, clustnum = 1, useCentroidDist = TRUE, vcmode = TRUE, outDir="./outputs"){
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # do ward's clustering on spatial and non-spatial attributes.
  # input:
  # (1) adata: attribute data (data.frame)
  # (2) pdata: polygons (list), cooridnates must be projected (which is CONST_projected_proj4string, e.g. x, y, rather than log and lat)
  # (3) ianmwh: interested attribute data column names and weights (data.frame) e.g.d = data.frame(cbind(ATTR_NAME=c("Morans_I","attr1","attr2"), ATTR_WEIGHT=c(0.4,0.3,0.3)))
  # (4) snswh: spatial and non-spatial weights (vector)
  # (5) dthresh: geo distance threshold (numeric), kilometers. If the distance between two polygons exceeds this value, they cannot be merged into one cluster
  # (6) proj4string: output shp file projection information, so the output can be reprojected or no-projection applied (using GCS)
  # (7) clustnum: target clustering number
  # (8) useCentroidDist: if true, use the centroid distance instead of the hausdorff distance to measure the distance between polygons. It speeds up the process enormously. 
  # (9) vcmode: value chain mode, if false, the value (job numbers) of interested attributes (job categories) will add up (weights NOT applied) as a new column ('vcvalue') and clustering will be performed on this column. If true, clustering will be performed on all interested attributes.
  # (10) outDir: the relative dir path for ouput data
  # output:
  # new shp file with a new data column marks its ward's cluster number
  #
  # algorithm process: to avoid duplicated pdf distance computation between unchanged polgyons, a pdf distance matrix (which is faster than dataframe) is initilized and 
  # maintained in each loop
  #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  print(paste("=== algorithm starts at ", Sys.time(), " ==="))
  # do some inputs validation
  
  # start time
  algStartTime = Sys.time()
  
  # assign Interested Non-Spatial Attribute Names
  IN_NS_ATTR_NAMES = as.character(ianmwh[[1]])
  
  # assign Interested Non-Spatial Attribute Weights
  IN_NS_ATTR_WEIGHTS = as.numeric(as.character(ianmwh[[2]]))
  
  # assign distance threshold
  CONST_DISTANCE_THRESHOLD = dthresh
  
  # append a column to hold the original polygon/data sequence
  adata[,"orgIdx"] = c(1:nrow(adata))
  
  # append a column to hold the ward's cluster number
  adata[,"wardclut"] = c(1:nrow(adata))
  
  # append a column to hold the add-up value
  if(vcmode==FALSE){
    if(length(IN_NS_ATTR_NAMES)>1){
      adata[,"vcvalue"] = rowSums(adata[,IN_NS_ATTR_NAMES])
    } else {
      adata[,"vcvalue"] = adata[,IN_NS_ATTR_NAMES]
    }
  }
  
  # backup data, adata and pdata will be changed during clustering process
  adata_bak = adata
  pdata_bak = pdata
  
  # do the normalization on each attribute column before calculating the distance, NA will be retained in the result
  norm_ns_attrs_dataframe = adata[,IN_NS_ATTR_NAMES, drop=FALSE]
  for(ns_attr_no in 1:length(IN_NS_ATTR_NAMES)){
    norm_ns_attrs_dataframe[[ns_attr_no]] = f_norm(norm_ns_attrs_dataframe[[ns_attr_no]])
  }
  
  if(vcmode==FALSE){
    # do the normalization on "vcvalue" column
    norm_ns_attrs_dataframe[,"vcvalue"] = f_norm(adata[,"vcvalue"])
  }
  
  # get the current polygon numbers
  CUR_POLYGON_NUM = nrow(adata)
  
  # init a matrix (dist_mat) to record non-normalized s_dist and ns_dist, as well as the polygon idx pair: idx_i, idx_j
  # in each merge step, the len_dist_mat will be updated 
  # at the very beginning, the matrix contains n*(n-1)/2 elements with default values (all set to 0), 
  
  # set the row number.  
  len_dist_dataframe = CUR_POLYGON_NUM * (CUR_POLYGON_NUM - 1) / 2
  # init distance (spatial and non-spatial) dataframe with all 0 values
  # it is much faster to pre-allocate memory for the matrix and then update its elements than to use "rbind" to append new rows to an existing matrix 
  dist_dataframe = data.frame(cbind(s_dist=rep(0,len_dist_dataframe), ns_dist=rep(0,len_dist_dataframe), idx_i=rep(0,len_dist_dataframe), idx_j=rep(0,len_dist_dataframe)))
  # convert into matrix to boost processing speed. 
  dist_mat = as.matrix(dist_dataframe) 
  dist_mat_rowcounter = 1
  
  debugPrint(paste("initialization distance matrix starts at ", Sys.time()))
  
  # assign values to distance matrix
  for (idx_i in 1:(CUR_POLYGON_NUM-1)){
    latlon_i = pdata[[idx_i]]@labpt
    pi = SpatialPolygons(list(pdata[[idx_i]]))
    len_pi = gLength(pi)
    
    debugPrint(sprintf("handling p%d of %d",idx_i, CUR_POLYGON_NUM))
    
    for (idx_j in (idx_i+1):CUR_POLYGON_NUM){
      latlon_j = pdata[[idx_j]]@labpt
      
      if (useCentroidDist == TRUE){
        sp_dist = f_getCentroidDistance(latlon_i, latlon_j)
      } else {
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
        #debugPrint(sprintf("(p%d,p%d) hdist_adjust: %.2f,    shared: %.2f,    hdist: %.2f",idx_i, idx_j, sp_dist, len_shared_boundary_pij, hdist))
      }
      
      # test on the geo distance threshold, ignore if exceeds
      if (sp_dist > CONST_DISTANCE_THRESHOLD * 1000) next
      
      # calc non spatial distance (Euclidian distance, or Manhattan distance)
      # pre-conditions: (1) attribute is numeric (2) attribute is normalized (i.e., 0-1)
      # issue: how to handle the distance between NA attribute value(s)? A simple way is to set the distance 0 so it makes no contribution (not exactly, it actually shortens the "real" distance) to the final dPDF.
      # updates: if NA exists in any one of the calculating attributes, the resulst is NA, and will be retained in the raw dataframe
      nsp_dist = 0
      if(vcmode==FALSE){
        nsp_dist = abs(adata[idx_i, "vcvalue"] - adata[idx_j, "vcvalue"])
      } else {
        for(ns_attr_no in 1:length(IN_NS_ATTR_NAMES)){
          nsp_dist = nsp_dist + IN_NS_ATTR_WEIGHTS[ns_attr_no]*(norm_ns_attrs_dataframe[idx_i, ns_attr_no] - norm_ns_attrs_dataframe[idx_j, ns_attr_no])^2
        }
        nsp_dist = sqrt(nsp_dist)
      }
      
      
      # update row values for dist_mat
      dist_mat[dist_mat_rowcounter,"s_dist"] = sp_dist
      dist_mat[dist_mat_rowcounter,"ns_dist"] = nsp_dist
      dist_mat[dist_mat_rowcounter,"idx_i"] = idx_i
      dist_mat[dist_mat_rowcounter,"idx_j"] = idx_j
      dist_mat_rowcounter = dist_mat_rowcounter + 1
      
    }
  }
  
  # there may be rows in the matrix not getting updated because of the CONST_DISTANCE_THRESHOLD condition, remove them
  filter = dist_mat[,"idx_i"] > 0
  dist_mat = dist_mat[filter, , drop=FALSE]
  
  debugPrint(sprintf("initial dist_mat row number: %i",nrow(dist_mat)))
  debugPrint(paste("initialization distance matrix ends at ", Sys.time()))
  
  # perform ward's clustering
  while(nrow(adata) > 1) #continue if more than one rows exist
  {
    # reset min_idx_i and min_idx_j. They are always the index of orignal adata (i.e. adata_bak)
    min_idx_i = -1
    min_idx_j = -1
    
    # get the current polygon numbers
    CUR_POLYGON_NUM = nrow(adata)
    
    # if target cluster number reached, stop algorithm
    if (CUR_POLYGON_NUM == clustnum){
      print("=== Target cluster number reached, exit. ===")
      break
    }
    
    debugPrint(sprintf("merging, %d polygons remain",CUR_POLYGON_NUM))
    
    # make a workable copy of dist_mat in each loop
    filtered_dist_mat = dist_mat
        
    if (nrow(filtered_dist_mat) == 0){ # means the geodistance threshold filter out all polggon pairs, it's time to exit
      print("=== No more polgyons can be merged, exit. ===")
      break
    }
    
    # normalize the s_dist (ns_dist) before computing the dPDF
    filtered_dist_mat[,"s_dist"] = f_norm(filtered_dist_mat[,"s_dist"])
    # the meaning of using normalized distance value on 'normalized attributes' is only for distance value comparison 
    # NA is set to 0.5 as the final normalized distance value
    filtered_dist_mat[,"ns_dist"] = f_norm(filtered_dist_mat[,"ns_dist"], na.rep.auto = FALSE, na.rep = 0.5)
    
    
    # to find the smallest dPDF in filtered_dist_mat, put the dPDF result in "s_dist" to save some memory
    filtered_dist_mat[,"s_dist"] = snswh[1]*filtered_dist_mat[, "s_dist"] + snswh[2]*filtered_dist_mat[, "ns_dist"]
    # sort filtered_dist_mat on "s_dist" column
    filtered_dist_mat = filtered_dist_mat[order(filtered_dist_mat[,"s_dist"]), ,drop=FALSE]
    # get the original polygon index wantted
    min_idx_i = min(as.integer(filtered_dist_mat[1,"idx_i"]), as.integer(filtered_dist_mat[1,"idx_j"]))
    min_idx_j = max(as.integer(filtered_dist_mat[1,"idx_i"]), as.integer(filtered_dist_mat[1,"idx_j"]))
    min_pdf_dist = filtered_dist_mat[1,"s_dist"]
    
    debugPrint(sprintf("min pdf distance found between (p%d,p%d) :  %.10f",min_idx_i, min_idx_j, min_pdf_dist))
    
    idx_ij_filter = (adata[["orgIdx"]] == min_idx_i) | (adata[["orgIdx"]] == min_idx_j)
    idx_i_filter = adata[["orgIdx"]] == min_idx_i
    idx_j_filter = adata[["orgIdx"]] == min_idx_j
    
    # merge two polygons: p(min_idx_i) and p(min_idx_j)
    # step1. handle the geometry data
    tmp_pij = SpatialPolygons(pdata[idx_ij_filter])
    uni_pij = unionSpatialPolygons(tmp_pij, c(1,1))
    uni_pij@polygons[[1]]@ID = paste(GLOBAL_polygon_id_prefix, as.character(GLOBAL_polygon_id_counter), sep="")
    GLOBAL_polygon_id_counter <<- GLOBAL_polygon_id_counter + 1 #change variable outside the function scope
    
    # store the newly merged polygon to Polygon_i
    pdata[idx_i_filter] = uni_pij@polygons
        
    # step2. handle the attribute data
    # merge attribute data into min_idx_i, how to handle NA should be considered
    ns_attr_no = 1
    for(ns_attr_no in 1:length(IN_NS_ATTR_NAMES)){
      adata[idx_i_filter, IN_NS_ATTR_NAMES[ns_attr_no]] = adata[idx_i_filter, IN_NS_ATTR_NAMES[ns_attr_no]] +  adata[idx_j_filter, IN_NS_ATTR_NAMES[ns_attr_no]]
    }
    if(vcmode==FALSE){
      adata[idx_i_filter,"vcvalue"] = adata[idx_i_filter,"vcvalue"] + adata[idx_j_filter,"vcvalue"]
    }
    
    # find the wardclust value of row min_idx_j in dataTargetYear before remove it
    wardclust_j_val = adata[idx_j_filter, "wardclut"]
    # find the wardclust value that will be used to replace wardclust_j_val in the dataTargetYear_bak
    wardclust_i_val = adata[idx_i_filter, "wardclut"]
        
    # update wardclust column in the dataTargetYear_bak
    # find all rows in dataTargetYear_bak whose wardclust value == wardclust_j_val
    tmpfilter = adata_bak[, "wardclut"] == wardclust_j_val
    adata_bak[tmpfilter, "wardclut"] = wardclust_i_val
    
    # then remove Polygon_j from the list
    pdata[idx_j_filter] = NULL
    # remove adata row on min_idx_j
    adata = adata[!idx_j_filter,]
    # recalculate the idx_i_filter because the length of adata is changed
    idx_i_filter = adata[["orgIdx"]] == min_idx_i    
    
    # if there are only 2 ploygons, and when reach this line, they have been merged, ignore the update distance data frame process
    if (CUR_POLYGON_NUM == 2) next
    
    # update the distance mat by removing all rows whose idx_i is min_idx_i or min_idx_j, or whose idx_j is min_idx_i or min_idx_j
    dist_mat_filter = dist_mat[,"idx_i"]!=min_idx_i & dist_mat[,"idx_i"]!=min_idx_j & dist_mat[,"idx_j"]!=min_idx_i & dist_mat[,"idx_j"]!=min_idx_j
    dist_mat = dist_mat[dist_mat_filter, , drop=FALSE]

    norm_ns_attrs_dataframe = adata[,IN_NS_ATTR_NAMES, drop=FALSE]
    
    # then do the normalization, NA will be retained in the result
    for(ns_attr_no in 1:length(IN_NS_ATTR_NAMES)){
      norm_ns_attrs_dataframe[[ns_attr_no]] = f_norm(norm_ns_attrs_dataframe[[ns_attr_no]])
    }
    
    if(vcmode==FALSE){
      # do the normalization on "vcvalue" column
      norm_ns_attrs_dataframe[,"vcvalue"] = f_norm(adata[,"vcvalue"])
    }
    
    # newly created polygon
    latlon_n = pdata[idx_i_filter][[1]]@labpt
    pN = SpatialPolygons(pdata[idx_i_filter])
    len_pN = gLength(pN)
    
    # rest polygon//data list
    pRest = pdata[!idx_i_filter]
    aRest = adata[!idx_i_filter,]
    norm_ns_attrs_dataframe_rest = norm_ns_attrs_dataframe[!idx_i_filter,]
        
    tmp_len = length(pRest)
    tmp_dist_dataframe = data.frame(cbind(s_dist=rep(0,tmp_len), ns_dist=rep(0,tmp_len), idx_i=rep(0,tmp_len), idx_j=rep(0,tmp_len)))
    tmp_dist_mat = as.matrix(tmp_dist_dataframe) 
    tmp_rowcounter = 1
    
    # append new rows in dist_mat
    for (rest_idx in 1:length(pRest)){
      latlon_r = pRest[[rest_idx]]@labpt
      
      if (useCentroidDist == TRUE){
        sp_dist = f_getCentroidDistance(latlon_n, latlon_r)
      } else {  
        # calc spatial distance
        pR = SpatialPolygons(pRest[rest_idx])
        len_pR = gLength(pR)
        hdist =  gDistance(pN,pR,hausdorff=TRUE)
        
        pNR = SpatialPolygons(list(pdata[idx_i_filter][[1]],pRest[[rest_idx]]))
        union_pNR = unionSpatialPolygons(pNR ,c(1,1))
        len_union_pNR = gLength(union_pNR)
        len_shared_boundary_pNR = len_pN + len_pR - len_union_pNR
        
        # ref: Joshi' paper p. 22
        sp_dist = hdist * (1 - 2 * len_shared_boundary_pNR / (len_pN + len_pR))
      }
      
      if (sp_dist > CONST_DISTANCE_THRESHOLD * 1000) next
      
      # calc non spatial distance (Euclidian distance, or we can use Manhattan distance)
      nsp_dist = 0
      if(vcmode==FALSE){
        nsp_dist = abs(norm_ns_attrs_dataframe[idx_i_filter, "vcvalue"] - norm_ns_attrs_dataframe_rest[rest_idx, "vcvalue"])
      } else {
        for(ns_attr_no in 1:length(IN_NS_ATTR_NAMES)){
          nsp_dist = nsp_dist + IN_NS_ATTR_WEIGHTS[ns_attr_no]*(norm_ns_attrs_dataframe[idx_i_filter, ns_attr_no] - norm_ns_attrs_dataframe_rest[rest_idx, ns_attr_no])^2
        }
        nsp_dist = sqrt(nsp_dist)
      }
      # set row values for dist_dataframe
      tmp_dist_mat[tmp_rowcounter,"s_dist"] = sp_dist
      tmp_dist_mat[tmp_rowcounter,"ns_dist"] = nsp_dist
      tmp_dist_mat[tmp_rowcounter,"idx_i"] = min(min_idx_i, aRest[rest_idx,"orgIdx"])
      tmp_dist_mat[tmp_rowcounter,"idx_j"] = max(min_idx_i, aRest[rest_idx,"orgIdx"])
      tmp_rowcounter = tmp_rowcounter + 1
    }
    
    filter = tmp_dist_mat[,"idx_i"] > 0
    tmp_dist_mat = tmp_dist_mat[filter, , drop=FALSE]
    
    dist_mat = rbind(dist_mat, tmp_dist_mat)
}
  
  # renumbering the cluster ID in sequence
  orgclutIds = as.integer(levels(factor(adata[,"wardclut"])))
  orgclutIds = sort(orgclutIds)
  for(i in 1:length(orgclutIds)){
    filter = adata[,"wardclut"] == orgclutIds[i]
    adata[filter,"wardclut"] = i
    
    filter = adata_bak[,"wardclut"] == orgclutIds[i]
    adata_bak[filter,"wardclut"] = i
  }
  
  # save merged result into a tmp shp file
  sp = SpatialPolygons(pdata)
  sp@proj4string = CRS(CONST_projected_proj4string)
  # for the merged polygons, the attribute data should only contain the job numbers and cluster number as valid information
  if(vcmode==FALSE){
    adata = adata[,c(IN_NS_ATTR_NAMES,"wardclut","vcvalue")]
  } else {
    adata = adata[,c(IN_NS_ATTR_NAMES,"wardclut")]
  }
  newDataFrame = SpatialPolygonsDataFrame(sp,data=adata, match.ID = FALSE)
  newDataFrame_pj = spTransform(newDataFrame,CRS(CONST_EPSG4283_proj4string))
  
  # save a new copy of origial polygons with updated cluseter information
  sp_bak =SpatialPolygons(pdata_bak)
  sp_bak@proj4string = CRS(CONST_projected_proj4string)
  newDataFrame_bak = SpatialPolygonsDataFrame(sp_bak,data=adata_bak, match.ID = FALSE)
  newDataFrame_pj_bak = spTransform(newDataFrame_bak,CRS(CONST_EPSG4283_proj4string))
  
  #plot(newDataFrame_pj_bak)
  #plot(newDataFrame_pj)
  
  # output data in shapefile format
  writeOGR(obj=newDataFrame_pj, dsn=outDir, layer="dump_clustered", driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
  writeOGR(obj=newDataFrame_pj_bak, dsn=outDir, layer="dump_original", driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
  
  # output data(non geometry part) in stata format
  stataFilePath = sprintf("%s/dump_clustered.dta", outDir)
  statabakFilePath = sprintf("%s/dump_original.dta", outDir)
  if(file.exists(stataFilePath)) file.remove(stataFilePath)
  if(file.exists(statabakFilePath)) file.remove(statabakFilePath)
  write.dta(newDataFrame_pj@data, stataFilePath)
  write.dta(newDataFrame_pj_bak@data, statabakFilePath)
  
  # output a png file for quick review
  f_visualize_dataframe(newDataFrame_pj_bak, pngFilePath=sprintf("%s.png", outDir), geodisthreshold=dthresh, spatialWeight=snswh[1])
  
  # output data in geojson format
  #geoJSONFilePath = sprintf("%s/dump_clustered.geojson", outDir)
  #geoJSONbakFilePath = sprintf("%s/dump_original.geojson", outDir)
  # for spatial data frame, convert to geojson
  #if(file.exists(geoJSONFilePath)) file.remove(geoJSONFilePath)
  #if(file.exists(geoJSONbakFilePath)) file.remove(geoJSONbakFilePath)
  #writeOGR(obj=newDataFrame_pj, dsn=geoJSONFilePath, layer="tmp", driver = "GeoJSON",  check_exists=TRUE, overwrite_layer=TRUE)    
  #writeOGR(obj=newDataFrame_pj_bak, dsn=geoJSONbakFilePath, layer="tmp", driver = "GeoJSON",  check_exists=TRUE, overwrite_layer=TRUE)
  
  # output data in KML format
  #kmlFilePath = sprintf("%s/dump_clustered.kml", outDir)
  #kmlbakFilePath = sprintf("%s/dump_original.kml", outDir)
  # for spatial data frame, convert to geojson
  #if(file.exists(kmlFilePath)) file.remove(kmlFilePath)
  #if(file.exists(kmlbakFilePath)) file.remove(kmlbakFilePath)
  #writeOGR(obj=newDataFrame_pj, dsn=kmlFilePath, layer="tmp", driver = "KML",  check_exists=TRUE, overwrite_layer=TRUE)    
  #writeOGR(obj=newDataFrame_pj_bak, dsn=kmlbakFilePath, layer="tmp", driver = "KML",  check_exists=TRUE, overwrite_layer=TRUE)
  
  algEndTime = Sys.time()
  print(paste("=== algorithm ends at ", Sys.time(), " ==="))
  print(sprintf("=== all done (in %.2f seconds) ===", as.numeric(algEndTime-algStartTime, units="secs")))
}

f_run <- function(num = -1, 
                  geodisthreshold = 10, 
                  targetclusternum = 1, 
                  useCentroidDist = TRUE,
                  displayColNames = c("LGA_CODE", "LGA", "ZONE_CODE"),
                  interestedColNames = c("X2310"),
                  interestedColWeights = c(1),
                  spatialNonSpatialDistWeights = c(0.5, 0.5),
                  ignoreEmptyRow = TRUE,
                  ignoreEmptyRowJobNum = 1,
                  vcmode = TRUE, 
                  outDir = "./outputs"
                  ){
  # input:
  # (1) num: number of rows taken from the original dataset, -1 means the entire data rows will be used
  # (2) geodisthreshold: geo distance threshold (numeric), kilometers. If the distance between two polygons exceeds this value, they cannot be merged into one cluster
  # (3) targetclusternum: target clustering number
  # (4) useCentroidDist: if true, use the centroid distance instead of the hausdorff distance to measure the distance between polygons. It speeds up the process enormously. 
  # (5) displayColNames: the columns will appear in the result shape file, no overlap with interestedColNames
  # (6) interestedColNames: interested attribute data column names, 
  # (7) interestedColWeights: interested attribute data column weights, should sum up to 1
  # (8) spatialNonSpatialDistWeights: spatial and non-spatial weights (vector)
  # (9) ignoreEmptyRow: if true, input dataset will be downsized by ignoreEmptyRowJobNum filter. 
  # (10) ignoreEmptyRowJobNum: column value less than this threshold is treated as empty row
  # (11) vcmode: value chain mode, if false, the value (job numbers) of interested attributes (job categories) will add up (weights NOT applied) as a new column ('vcvalue') and clustering will be performed on this column. If true, clustering will be performed on all interested attributes.
  # (12) outDir: the relative dir path for ouput data
  
  # use fast read mode
  #x <- readShapePoly("./data/SplitPoly_X_Employment_fullcode")
  x <- readShapePoly(sprintf("./data/MidPoly_X_Employment_%s",gDataYear))
  # assign proj infomation
  attr(x@proj4string,"projargs") = "+proj=longlat +ellps=GRS80 +no_defs"
  original_proj4string = attr(x@proj4string,"projargs")
  
  # uppercase all column names
  colnames(x@data) = toupper(colnames(x@data))
  
  # check if the original data is projected
  # Transform the polygons (which were read in as unprojected geographic coordinates) to an Albers Equal Area projection
  # this is essential since we need calculate distance/area for polygons
  if (is.projected(x)){
    x_pj = x
  } else {
    x_pj = spTransform(x,CRS(CONST_projected_proj4string))
  }
  
  if (num <= 0 | num > nrow(x_pj@data)){
    DATA_ROW_NUM = nrow(x_pj@data)
  } else {
    DATA_ROW_NUM = num
  }
  
  # uppercase all display and interested column names
  displayColNames = toupper(displayColNames)
  interestedColNames = toupper(interestedColNames)
  
  testData = x_pj@data[1:DATA_ROW_NUM, c(displayColNames,interestedColNames)]
  testData[,"wardclut"] = 1:DATA_ROW_NUM
  testPolyList = x_pj@polygons[1:DATA_ROW_NUM]
  
  # ignore rows if all interested column values are 0
  if (ignoreEmptyRow==TRUE){
    filter = rep(FALSE, nrow(testData))
    for(colname in interestedColNames){
      filter = filter | (testData[,colname] >= ignoreEmptyRowJobNum)
    }
    testData = testData[filter,]
    testPolyList = testPolyList[filter]
  }
  
  nmwt=data.frame(cbind(ATTR_NAME=interestedColNames, ATTR_WEIGHT=interestedColWeights))
  
  if (outDir == "./outputs"){
    outDir = paste(outDir,gDataYear,sep="/")
  }
  f_wards(adata=testData, pdata=testPolyList, ianmwh=nmwt, snswh=spatialNonSpatialDistWeights, dthresh=geodisthreshold, proj4string=original_proj4string, clustnum=targetclusternum, useCentroidDist=useCentroidDist, vcmode = vcmode, outDir=outDir)
}

f_visualize_dataframe <- function(map = NULL,
                        pngFilePath = "./tmp.png",
                        withLegend = FALSE,
                        geodisthreshold = 10,
                        spatialWeight = 0.1
) {
  # render clusters in colors on the original polygons
  if(is.null(map)) return(FALSE)

  clusterIds = as.numeric(levels(factor(map@data$wardclut)))
  png(file=pngFilePath, width=800, height=800)
  plot(map, col=gColVec[map@data$wardclut], axes=F)
  if(withLegend==TRUE){
    legend(x=map@bbox["x","min"],y=map@bbox["y","max"], legend=as.character(clusterIds), fill=gColVec[clusterIds])
  }
  title(sprintf("[%s] Ploygons: %i, Clusters: %i, Distance: %.2fKM, SpatialWeight: %.1f", gDataYear, nrow(map@data), length(clusterIds),geodisthreshold, spatialWeight))
  dev.off()
  return(TRUE)
}

# this function can be called separately to visualize the clusters
f_visualize <- function(shpFilePath="./outputs/dump_clustered", 
                        withLegend=FALSE,
                        geodisthreshold = 10,
                        spatialWeight = 0.1
                        ) {
  # render clusters in colors on the original polygons
  if(!file.exists(sprintf("%s.shp",shpFilePath))) return(FALSE)
  
  map <- readShapePoly(shpFilePath)
  clusterIds = as.numeric(levels(factor(map@data$wardclut)))
  png(file=sprintf("%s.png",shpFilePath), width=800, height=800)
  #par(mfrow=c(1,2))
  #plot(map, col=gColVec[map@data$wardclut], axes=F)
  plot(map, col=gColVec[map@data$wardclut], axes=F)
  if(withLegend==TRUE){
    legend(x=map@bbox["x","min"],y=map@bbox["y","max"], legend=as.character(clusterIds), fill=gColVec[clusterIds])
  }
  title(sprintf("Clusters: %i, Distance: %.2fKM, SpatialWeight: %.1f", length(clusterIds),geodisthreshold, spatialWeight))
  dev.off()
  return(TRUE)
}


f_run_globaldata <- function(num = -1, 
                  geodisthreshold = 10, 
                  targetclusternum = 1, 
                  useCentroidDist = TRUE,
                  displayColNames = c("LGA_CODE", "LGA", "ZONE_CODE"),
                  interestedColNames = c("X2310", "X2412", "X8500"),
                  interestedColWeights = c(0.333, 0.333, 0.333),
                  spatialNonSpatialDistWeights = c(0.5, 0.5),
                  ignoreEmptyRow = TRUE,
                  ignoreEmptyRowJobNum = 10,
                  vcmode = TRUE, 
                  outDir = "./outputs"
){
  # use fast read mode
  f_init_globaldata()
  
  if (num <= 0 | num > nrow(gX_pj@data)){
    DATA_ROW_NUM = nrow(gX_pj@data)
  } else {
    DATA_ROW_NUM = num
  }
  
  # uppercase all display and interested column names
  displayColNames = toupper(displayColNames)
  interestedColNames = toupper(interestedColNames)
  
  testData = gX_pj@data[1:DATA_ROW_NUM, c(displayColNames,interestedColNames)]
  testData[,"wardclut"] = 1:DATA_ROW_NUM
  testPolyList = gX_pj@polygons[1:DATA_ROW_NUM]
  
  # ignore rows if all interested column values are 0
  if (ignoreEmptyRow==TRUE){
    filter = rep(FALSE, nrow(testData))
    for(colname in interestedColNames){
      filter = filter | (testData[,colname] >= ignoreEmptyRowJobNum)
    }
    testData = testData[filter,]
    testPolyList = testPolyList[filter]
  }
  
  nmwt=data.frame(cbind(ATTR_NAME=interestedColNames, ATTR_WEIGHT=interestedColWeights))
  if (outDir == "./outputs"){
    outDir = paste(outDir,gDataYear,sep="/")
  }
  f_wards(adata=testData, pdata=testPolyList, ianmwh=nmwt, snswh=spatialNonSpatialDistWeights, dthresh=geodisthreshold, proj4string=gOriginal_proj4string, clustnum=targetclusternum, useCentroidDist=useCentroidDist, vcmode = vcmode, outDir=outDir)
}


f_batch <- function(vecGeoDist = c(1,5,10,20),
                    vecSpatialWeights = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                    interestedColNames = c("X2310"),
                    interestedColWeights = c(1),
                    ignoreEmptyRowJobNum = 1,
                    vcmode = TRUE,
                    outDir = "X2310"
                    ){
  
  for(geodist in vecGeoDist){
    for(sw in vecSpatialWeights){
      fullOutputDir = sprintf("./outputs/%s/%s/%.2f/%.1f", gDataYear, outDir, geodist, sw)
      print(sprintf("=== run for %s", fullOutputDir))
      dir.create(fullOutputDir, showWarnings=FALSE, recursive=TRUE)
      f_run_globaldata(interestedColNames=interestedColNames, 
            interestedColWeights=interestedColWeights,
            ignoreEmptyRowJobNum=ignoreEmptyRowJobNum,
            vcmode=vcmode,
            outDir=fullOutputDir,
            geodisthreshold = geodist,
            spatialNonSpatialDistWeights = c(sw,(1-sw))
            ) 
    }
  }
}

f_init_globaldata <-function(){
  if(is.null(gX)){
    #gX <<- readShapePoly("./data/SplitPoly_X_Employment_fullcode")
    gX <<- readShapePoly(sprintf("./data/MidPoly_X_Employment_%s",gDataYear))
    # assign proj infomation
    attr(gX@proj4string,"projargs") <<- "+proj=longlat +ellps=GRS80 +no_defs"
    gOriginal_proj4string <<- attr(gX@proj4string,"projargs")
    
    # uppercase all column names
    colnames(gX@data) = toupper(colnames(gX@data))
    
    # check if the original data is projected
    # Transform the polygons (which were read in as unprojected geographic coordinates) to an Albers Equal Area projection
    # this is essential since we need calculate distance/area for polygons
    if (is.projected(gX)){
      gX_pj <<- gX
    } else {
      gX_pj <<- spTransform(gX,CRS(CONST_projected_proj4string))
    }
  }
}
