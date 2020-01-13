extract_year <- function (coords, data_year, files, year_vec){




}


extract_year_month <- function (coords, data_year, data_month, files, year_vec, month_vec){

  cov <- rep(NA, nrow(coords))
  
  unq <- distinct(data.frame(data_year, data_month))
  
  ext <- function(u){
    # Year might not exists, so take nearest.
    if(unq$data_year[u] %in% year_vec){
      year <- unq$data_year[u]
    } else {
      year <- year_vec[which.min(abs(year_vec - unq$data_year[u] - 0.0001))]
    }
    
    # Find the right flie and read raster
    file_ii <- year_vec == year & month_vec == unq$data_month[u]
    stopifnot(sum(file_ii) == 1)
    r <- raster(files[file_ii])
    
    vals <- extract(r, coords[data_year == unq$data_year[u]  & data_month == unq$data_month[u], ])
    return(vals)
  }
  
  val_list <- mclapply(seq_len(nrow(unq)), ext, mc.cores = detectCores())
  
  # For each unique month year combination, read that raster and extract.
  for(u in length(val_list)){
    # Extract
    cov[data_year == unq$data_year[u] & data_month == unq$data_month[u]] <-
      val_list[[u]]
  }

  return(cov)


}