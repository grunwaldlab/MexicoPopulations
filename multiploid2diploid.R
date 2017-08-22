multiploid2diploid <- function(x, to = 2){
  
  ploidy_tab <- info_table(x, type = "ploidy")
  mat_sum <- rowSums(ploidy_tab > to, na.rm = TRUE)
  
  diploid_pinf <- x[mat_sum == 0]
  
  diploid_pinf <- recode_polyploids(diploid_pinf, newploidy = to)
  
  return(diploid_pinf)
  
}


