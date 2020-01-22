
get_preds <- function(t){
  stopifnot(inherits(t, 'train'))
  
  oos <- t$pred %>% arrange(rowIndex)
  row_matches <- sapply(1:length(t$bestTune), function(x) oos[, names(t$bestTune)[x]] == t$bestTune[[x]])
  best_rows <- rowMeans(row_matches) == 1
  
  d <- oos[best_rows, ]
  
  return(as.numeric(d$pred))
  
} 

