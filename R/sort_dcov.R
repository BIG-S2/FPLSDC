sort_dcov <- function(X,img){
  dcov_vec <- NULL
  X <- scale(X)
  n <- dim(X)[1]
  dim_len <- dim(X)[2]
  mar_dcov <- rep(-1,dim_len)
  if (dim_len > floor(3*n/log(n))){
  d_img <- dist(img)
  for (i in 1:length(X[1,])) {
    dcov_vec[i] <- dcov(dist(X[,i]),d_img)
  }
  dcov_vec_sort <- sort(dcov_vec,decreasing = TRUE, index.return = TRUE)

  dcov_vec_sort_sub <- list(x = dcov_vec_sort$x[1:floor(3*n/log(n))],ix = dcov_vec_sort$ix[1:floor(3*n/log(n))])
  for (i in 1:floor(3*n/log(n))) {
    m_dcov <- marginal_dcov(X[,dcov_vec_sort_sub$ix[i]],img)
    mar_dcov[dcov_vec_sort_sub$ix[i]] <- -m_dcov$img_par$min_obj
  }
  sort_d <- sort(mar_dcov,decreasing = TRUE,index.return=TRUE)

  sort_index <- c(sort_d$ix[1:floor(3*n/log(n))], dcov_vec_sort$ix[-seq(floor(3*n/log(n)))])
  }else{
    for (i in 1:dim_len) {
      m_dcov <- marginal_dcov(X[,i],img)
      mar_dcov[i] <- -m_dcov$img_par$min_obj
    }
    sort_d <- sort(mar_dcov,decreasing = TRUE,index.return=TRUE)
    sort_index <- sort_d$ix
  }
  return(rank = sort_index)
}
