library(numDeriv)
library(energy)
library(Rcpp)
library(nloptr)
library(MASS)

sin_matrix <- function(phi){
  Spherical <- diag(length(phi))
  for (i in seq(length(phi))) {
    Spherical[i,i] <- sin(phi[i])
  }
  return(Spherical)
}
cos_matrix <- function(phi){
  Spherical <- diag(length(phi))
  for (i in seq(length(phi))) {
    Spherical[i,i] <- cos(phi[i])
  }
  return(Spherical)
}
beta_vector <- function(phi){
  beta_vec <- NULL
  max_matrix <- sin_matrix(phi)
  for (i in seq(length(phi)+1)) {
    if(i == 1) beta_vec[1] <- cos(phi[1])
    if(i > 1 && i < (length(phi)+1)){
      sub_matrix <- max_matrix[1:(i-1),1:(i-1)]
      beta_vec[i] <- det(matrix(sub_matrix,i-1,i-1))*cos(phi[i])
    }
    if(i == (length(phi)+1)) beta_vec[length(phi) + 1] <- det(max_matrix)
  }
  return(beta_vec)
}
grad_beta_vector <- function(phi){
  grad_vec <- matrix(0,length(phi)+1,length(phi))
  #grad_vec <- diag(length(phi)+1)
  cos_m <- cos_matrix(phi)
  sin_m <- sin_matrix(phi)
  grad_vec[1,1] <- -sin_m[1]
  if(length(phi) == 1){
    grad_vec[1,1] <- -sin_m[1]
    grad_vec[2,1] <- cos_m[1]
  }
  if(length(phi)>1){
    for (j in seq(length(phi))) {
      for (i in seq(j,length(phi)+1)) {
        if(i == j && i > 1){
          grad_vec[i,j] <- -sin_m[i,i]*det(matrix(sin_m[seq(i-1),seq(i-1)],i-1,i-1))
        }else{
          if(i == 2){
            grad_vec[i,j] <- cos_m[i,i]*cos_m[j,j]
          }
          if(i > 2 && i < length(phi)+1){
            grad_vec[i,j] <- cos_m[i,i]*cos_m[j,j]*det(matrix(sin_m[seq(i-1)[-j],seq(i-1)[-j]],i-2,i-2))
          }
        }
        if(i == length(phi)+1){
          grad_vec[i,j] <- cos_m[j,j]*det(matrix(sin_m[seq(i-1)[-j],seq(i-1)[-j]],i-2,i-2))
        }
      }
    }
  }
  return(grad_vec)
}


scalar_cons <- function(x){
  return(t(x) %*% x - 1)
}
scalar_metrospolis <- function(f,f_new,T){
  if(abs(f - f_new) < 10^-6) return(0)
  if(f > f_new){
    return(1)
  }else{
    p <- exp((f - f_new)/T)
    if(p > runif(1)) return(1)
    else return(0)
  }
}
# the first opt
scalar_obj_sqp <- function(X,Y,index,beta){
  return(-dcovU(X[,index] %*% beta, Y))
}
scalar_grad1 <- function(X_sub, Y_sub, index, beta, n){
  grad_dcovU <- -grad_dcov_U(sign(B_distance_2(X_sub[,index] %*% beta_vector(beta))), B_distance(Y_sub), X_sub[,index] %*% grad_beta_vector(beta),rowSumsC(B_distance(Y_sub)),sum(rowSumsC(B_distance(Y_sub))),min(n,2000))
  return(grad_dcovU)
}
scalar_generate1 <- function(X_sub, Y_sub, index, beta1, n){
  beta_new <- beta1 - runif(1,0,pi) * scalar_grad1(X_sub, Y_sub, index, beta1, n)
  return(beta_new)
}

scalar_generate1_rand <- function(X_sub, Y_sub, index, beta1, beta2, n){
  rand <- runif(length(beta1),0,pi)
  beta_new <- beta1 -  runif(1,0,pi) * (scalar_grad1(X_sub, Y_sub, index, beta1, n) + sign(scalar_grad1(X_sub, Y_sub, index, beta2, n))*rand)
  return(beta_new)
}
scalar_run1 <- function(times,X,Y,X_sub, Y_sub,index,beta){
  n <- dim(X)[1]
  f_list <- NULL
  beta_count <- NULL
  counts <- 2
  #beta <- v_norm(rnorm(8))
  f_list[1] <- scalar_obj_sqp(X,Y,index,beta_vector(beta))
  beta_count[[1]] <- beta
  breaknum <- 0
  beta_new <- scalar_generate1(X_sub, Y_sub, index, beta, n)
  beta1 <- beta
  T = 10
  for (i in 1:times) {
    f <- scalar_obj_sqp(X,Y,index,beta_vector(beta))
    f_new <- scalar_obj_sqp(X,Y,index,beta_vector(beta_new))
    m_num <- scalar_metrospolis(f = f,f_new = f_new,T)
    if(m_num == 0){
      breaknum <- breaknum + 1
      if(breaknum > 5) {
        break
      }
      beta_new <- scalar_generate1_rand(X_sub, Y_sub, index, beta, beta1, n)
    }
    if(m_num == 1){
      beta_count[[counts]] <- beta_new
      f_list[counts] <- f_new
      beta1 <- beta
      beta <- beta_new
      counts <- counts + 1
      breaknum <- 1
      if(f > f_new){
        beta_new <- scalar_generate1(X_sub, Y_sub, index, beta, n)
      }else{
        beta_new <- scalar_generate1_rand(X_sub, Y_sub, index, beta, beta1, n)
      }

    }
    T=T*0.95
  }
  min_index <- which.min(f_list)
  return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
}

# the second opt
scalar_obj_sqp_newbeta <- function(X,Y,new_index,beta){
  return(-dcovU(X[,new_index] %*% beta, Y))
}
scalar_grad2 <- function(X_sub, Y_sub, new_index, beta, n){
  grad_dcovU <- -grad_dcov_U(sign(B_distance_2(X_sub[,new_index] %*% beta_vector(beta))), B_distance(Y_sub), X_sub[,new_index] %*% grad_beta_vector(beta),rowSumsC(B_distance(Y_sub)),sum(rowSumsC(B_distance(Y_sub))),min(n,2000))
  return(grad_dcovU)
}
scalar_generate2 <- function(X_sub, Y_sub, new_index, beta1, n){
  beta_new <- beta1 - runif(1,0,pi) * scalar_grad2(X_sub, Y_sub, new_index, beta1, n)
  return(beta_new)
}
scalar_generate2_rand <- function(X_sub, Y_sub, new_index, beta1,beta2,n){
  rand <- runif(length(beta1),0,pi)
  beta_new <- beta1 -  runif(1,0,pi) * (scalar_grad2(X_sub, Y_sub, new_index, beta1, n) + sign(scalar_grad2(X_sub, Y_sub, new_index, beta2, n))*rand)
  return(beta_new)
}
scalar_run2 <- function(times,X,Y,X_sub, Y_sub,new_index,beta,n){
  f_list <- NULL
  beta_count <- NULL
  counts <- 2
  f_list[1] <- scalar_obj_sqp_newbeta(X,Y,new_index,beta_vector(beta))
  beta_count[[1]] <- beta
  breaknum <- 0
  beta1 <- beta
  T=10
  beta_new <- scalar_generate2(X_sub, Y_sub, new_index, beta, n)
  for (i in 1:times) {
    f <- scalar_obj_sqp_newbeta(X,Y,new_index,beta_vector(beta))
    f_new <- scalar_obj_sqp_newbeta(X,Y,new_index,beta_vector(beta_new))
    m_num <- scalar_metrospolis(f = f,f_new = f_new,T)
    if(m_num == 0){
      breaknum <- breaknum + 1
      if(breaknum > 5) {
        break
      }
      beta_new <- scalar_generate2_rand(X_sub, Y_sub, new_index, beta,beta1, n)
    }
    if(m_num == 1){
      beta_count[[counts]] <- beta_new
      f_list[counts] <- f_new
      beta1 <- beta
      beta <- beta_new
      counts <- counts + 1
      breaknum <- 1
      if(f > f_new){
        beta_new <- scalar_generate2(X_sub, Y_sub, new_index, beta, n)
      }else{
        beta_new <- scalar_generate2_rand(X_sub, Y_sub, new_index, beta, beta1, n)
      }
    }
    T=T*0.9
  }
  min_index <- which.min(f_list)
  return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
}

######## function
norm_1 <- function(pro_par){
  return(pro_par/sqrt(as.numeric(t(pro_par)%*% pro_par)))
}
scalar_obj_BB <-function(X,Y,pro_par,theta,sort_index,select_num,add_i,index){
  return(-dcovU(cos(theta)*(X[,index] %*% norm_1(pro_par))+sin(theta)*X[,sort_index$ix[-c(1:select_num)][add_i]],Y))
}
scalar_obj_grad_BB <- function(X,Y,index,pro_par,sort_index,add_i,select_num,theta){
  n <- dim(X)[1]
  if (length(X[,1]) > 2000){
    sample_index <- sample(length(X[,1]),2000)
    grad_ratio <- -grad_dcovU_ratio(sign(B_distance_2(cos(theta)*(X[sample_index,index] %*% norm_1(pro_par))+sin(theta)*X[sample_index,sort_index$ix[-c(1:select_num)][add_i]])),
                                    B_distance_2(-sin(theta)*(X[sample_index,index] %*% norm_1(pro_par))+cos(theta)*X[sample_index,sort_index$ix[-c(1:select_num)][add_i]]),
                                    B_distance((Y)[sample_index,]),rowSumsC(B_distance((Y)[sample_index,])),sum(rowSumsC(B_distance((Y)[sample_index,]))),length(sample_index))
  }
  else{
    grad_ratio <- -grad_dcovU_ratio(sign(B_distance_2(cos(theta)*(X[,index] %*% norm_1(pro_par))+sin(theta)*X[,sort_index$ix[-c(1:select_num)][add_i]])),
                                    B_distance_2(-sin(theta)*(X[,index] %*% norm_1(pro_par))+cos(theta)*X[,sort_index$ix[-c(1:select_num)][add_i]]),
                                    B_distance(Y),rowSumsC(B_distance(Y)),sum(rowSumsC(B_distance(Y))),n)
  }
  return(grad_ratio)
}
scalar_grad_norm <- function(x){
  return(sqrt(t(x) %*% x))
}
############
# scalar_generate3_rand <- function(X,Y,index,pro_par,sort_index,add_i,select_num,beta1,beta2){
#   rand <- runif(length(beta1),0,pi)
#   beta_new <- beta1 - runif(1,0,pi) * (scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,beta1) + sign(scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,beta2))*rand)
#   return(beta_new)
# }
# scalar_generate3 <- function(X,Y,index,pro_par,sort_index,add_i,select_num,beta1){
#   beta_new <- beta1 - runif(1,0,pi) * scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,beta1)
#   return(beta_new)
# }
# scalar_run3 <- function(X,Y,times,pro_par,beta,sort_index,select_num,add_i){
#   f_list <- NULL
#   beta_count <- NULL
#   counts <- 2
#   f_list[1] <- scalar_obj_BB(X,Y,pro_par,beta,sort_index,select_num,add_i)
#   beta_count[[1]] <- beta
#   breaknum <- 0
#   beta_new <- scalar_generate3(X,Y,index,pro_par,sort_index,add_i,select_num,beta)
#   T = 1
#   for (i in 1:times) {
#     f <- scalar_obj_BB(X,Y,index,pro_par,beta,sort_index,add_i,select_num,add_i)
#     f_new <- scalar_obj_BB(X,Y,index,pro_par,beta_new,sort_index,add_i,select_num,add_i)
#     m_num <- scalar_metrospolis(f = f,f_new = f_new,T)
#     if(m_num == 0){
#       breaknum <- breaknum + 1
#       if(breaknum > 5) {
#         break
#       }
#       beta_new <- scalar_generate3_rand(X,Y,index,pro_par,sort_index,add_i,select_num,beta,beta1)
#     }
#     if(m_num == 1){
#       beta_count[[counts]] <- beta_new
#       f_list[counts] <- f_new
#       beta1 <- beta
#       beta <- beta_new
#       counts <- counts + 1
#       breaknum <- 1
#       if(f > f_new){
#         beta_new <- scalar_generate3(X,Y,index,pro_par,sort_index,add_i,select_num,beta)
#       }else{
#         beta_new <- scalar_generate3_rand(X,Y,index,pro_par,sort_index,add_i,select_num,beta,beta1)
#       }
#     }
#     T=T*0.9
#   }
#   if(counts == 1) {
#     f_list <- f
#     beta_count[[counts]] <- beta
#   }
#   min_index <- which.min(f_list)
#   return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
# }
###################
scalar_BB_opt <- function(X,Y,index,pro_par,sort_index,add_i,select_num,x0){
  alpha <- 0.025
  M = 5
  c1 <- 0.5
  beta <- 0.5
  epslion <- 1e-4
  k=1
  max_k <- 10
  store_x <- numeric()
  store_value <- numeric()
  store_x[1] <- x0
  store_value[1] <- scalar_obj_BB(X,Y,pro_par,store_x[1],sort_index,select_num,add_i,index)
  while (TRUE) {
    for (i in 1:10) {
      beta_3 <- store_x[k] - alpha*scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,store_x[k])
      if(scalar_obj_BB(X,Y,pro_par,beta_3,sort_index,select_num,add_i,index) < max(store_value[(k + 1 - min(k,M)):k]) - c1*alpha*scalar_grad_norm(scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,store_x[k]))^2){
        break
      }
      alpha <- beta * alpha
    }
    store_x[k+1] <- store_x[k] - alpha*t(scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,store_x[k]))
    s <- store_x[k+1] - store_x[k]
    y <- scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,store_x[k+1]) - scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,store_x[k])
    alpha <- as.numeric((t(s) %*% y)/(t(y) %*% y))
    if(is.null(alpha)) alpha <- 0.01
    if(alpha < 0.01) {alpha <- 0.01}
    if(alpha > 0.5) {alpha <- 0.5}
    k <- k+1
    if(k >= max_k) break
    store_value[k] <- scalar_obj_BB(X,Y,pro_par,store_x[k],sort_index,select_num,add_i,index)
    if(scalar_grad_norm(alpha * scalar_obj_grad_BB(X,Y,index,pro_par,sort_index,add_i,select_num,store_x[k])) < epslion) break
    if(abs(store_value[k] - store_value[(k-1)]) < 1e-5) break
  }
  return(list(par = store_x, value = -store_value))
}
############
# scalar_run4 <- function(times,X,Y,pro_par,beta,sort_index,select_num,add_i){
#   f_list <- NULL
#   #T <- 100
#   beta_count <- NULL
#   counts <- 2
#   f_list[1] <- scalar_obj_BB(X,Y,pro_par,beta,sort_index,select_num,add_i)
#   beta_count[[1]] <- beta
#   breaknum <- 0
#   beta_new <- beta
#   beta1 <- beta
#   T = 1
#   for (i in 1:times) {
#     f <- scalar_obj_BB(X,Y,pro_par,beta,sort_index,select_num,add_i)
#     opt_f <- scalar_BB_opt(X,Y,pro_par,sort_index,select_num,add_i,beta_new)
#     f_new <- opt_f$value[which.max(opt_f$value)]
#     m_num <- scalar_metrospolis(f = f,f_new = f_new,T)
#     if(m_num == 0){
#       breaknum <- breaknum + 1
#       if(breaknum > 5) {
#         break
#       }
#       beta_new <- scalar_generate3_rand(X,Y,index,pro_par,sort_index,add_i,select_num,beta,beta1)
#     }
#     if(m_num == 1){
#       beta_count[[counts]] <- opt_f$par[which.max(opt_f$value)]
#       f_list[counts] <- f_new
#       beta1 <- beta_new
#       beta <- opt_f$par[which.max(opt_f$value)]
#       beta_new <- scalar_generate3_rand(X,Y,index,pro_par,sort_index,add_i,select_num,beta,beta1)
#       counts <- counts + 1
#       breaknum <- 1
#     }
#     T=T*0.9
#     if(T < 0.01){
#       break
#     }
#   }
#   min_index <- which.max(f_list)
#   return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
# }
#######################
scalar_select <- function(X,Y,select_num,cut_feature,min_feature=5,sort_d){
  n <- length(as.matrix(X)[,1])
  feature_num <- length(as.matrix(X)[1,])
  sort_index <- sort_d
  feature_num_sub <- floor(min(cut_feature,length(as.matrix(X)[1,])))
  sort_index_sub <- list(x = sort_d$x[1:feature_num_sub],ix = sort_d$ix[1:feature_num_sub])
  index <- sort_index_sub$ix[1:select_num]
  cat("the initial index:",index,"\n")
  ###### Determines the initial parameter values
  sqp_par <- matrix(runif(select_num-1,0,pi))
  # cat("initial par:",beta_vector(sqp_par),"\n")
  # cat("initial dcov:",dcovU(X[,index]%*% beta_vector(sqp_par), Y),"\n")
  if (length(X[,1]) > 2000){
    sample_n <- sample(length(X[,1]),2000)
    X_sub <- X[sample_n,]
    Y_sub <- Y[sample_n,]
  }else{
    X_sub <- X
    Y_sub <- Y
  }
  new_index <- index
  new_beta <- sqp_par

  run1_opt <- scalar_run1(40,X,Y,X_sub, Y_sub,index,sqp_par)
  sqp_par_new <- run1_opt$min_beta
  # cat("optimize dcov:",dcovU(X[,index]%*% beta_vector(sqp_par_new), Y),"\n")
  sqp_par <- sqp_par_new
  ###### store the variation of dcovU of each variable
  index_variation <- rep(-1,feature_num)
  feature_num1 <- feature_num
  iterations <- 1
  end_dcov <- NULL  # store the iterated dcov
  end_dcov[1] <- 0
  while (TRUE) {
    index_variation <- rep(-1,feature_num1)
    d_cov <- dcovU(X[,index]%*%beta_vector(sqp_par),Y)
    sqp_par1 <- beta_vector(sqp_par)
    for (min_i in 1:length(index)) {
      pro_par <- sqp_par1[-min_i]/sqrt(sum((sqp_par1[-min_i])^2))
      index_variation[index[min_i]] <-  abs(d_cov - dcovU(as.matrix(X[,index[-min_i]])%*%pro_par,(Y)))
    }
    pro_par <- sqp_par1 / sqp_par1[1]
    add_length <- max((select_num+1),feature_num_sub)
    for (add_i in 1:length(seq(add_length-select_num))) {
      # cat(add_i,"\n")
      index_variation[sort_index$ix[-c(1:select_num)][add_i]] <- abs(tail(scalar_BB_opt(X,Y,index,pro_par,sort_index,add_i,select_num,0)$value,1) - d_cov)
      # cat("the selected num is:",select_num,"the number of the insert variable",add_i)
    }
    sort_index_sub <- sort(index_variation,decreasing = TRUE,index.return=TRUE)
    sort_index_old <- sort_index
    sort_index$ix[1:feature_num_sub] <- sort_index_sub$ix[1:feature_num_sub]


    # sort_index <- sort(index_variation,decreasing = TRUE,index.return=TRUE)
    new_index <- sort_index$ix[1:select_num]

    cat("the new updated index:",new_index,"\n")
    #cat("After the update, the position of the active variable??", match(feature_index,sort_index$ix),"\n")

    ##### After the new index is determined, Update the coefficient of new vaiable
    if (n > 2000){
      sample_n <- sample(n,2000)
      X_sub <- X[sample_n,]
      Y_sub <- Y[sample_n,]

    }else {
      X_sub <- X
      Y_sub <- Y
    }
    ini_beta <- matrix(v_norm(runif(select_num-1,0,pi)))
    # cat("initial dcov:",dcovU(X[,new_index]%*% beta_vector(ini_beta),Y),"\n")
    new_beta <- scalar_run2(40,X,Y,X_sub, Y_sub,new_index,ini_beta,n)$min_beta
    # cat("new dcov:",dcovU(X[,new_index]%*%beta_vector(new_beta),Y),"\n")
    end_dcov[iterations+1] <- dcovU(X[,new_index] %*%beta_vector(new_beta), Y)
    if (iterations > 10){
      cat("Reaches the maximum number of iterations\n")
      #if(true_index[1] != -1)return(list(index= new_index, position = sort_index$ix, value = dcovU(X[,new_index] %*% new_beta, Y_sub)))
      scalar_list <- list(par= beta_vector(new_beta),index= sort_index$ix[1:select_num],position = sort_index$ix,value = dcovU(X[,new_index] %*% beta_vector(new_beta) , Y))
      return(scalar_list)
      break
    }
    if ((end_dcov[iterations]/end_dcov[iterations+1] - 1)>0.01){
      cat("The dcovU decrease, stop!\n")
      return(list(par= beta_vector(sqp_par),index= sort_index_old$ix[1:select_num],position = sort_index_old$ix,value = dcovU(X[,index] %*% beta_vector(sqp_par), Y)))
      break
    }
    if (sum(sort(index) == sort(new_index)) == select_num){
      cat("The iteration stop\n")
      scalar_list <- list(par= beta_vector(new_beta),index= sort_index$ix[1:select_num],position = sort_index$ix,value = dcovU(X[,new_index] %*% beta_vector(new_beta), Y))
      return(scalar_list)
      break
    }
    index <- new_index
    iterations <- iterations + 1
    sqp_par <- new_beta
    feature_num_sub <- floor(feature_num_sub/min_feature)
  }
}

scalar_f_opt <- function(X,Y,i,sort_d,cut_feature){
  store_value <- scalar_select(X,Y,select_num = i,min_feature = 5,cut_feature = cut_feature,sort_d = sort_d)
  cat("the selected variable:",store_value$index,"\n")
  return(list(index = store_value$index, position = store_value$position,dcov = store_value$value,X_par = store_value$par))
}

scalar_fpls_select <- function(X,Y){
  Y <- matrix(Y,length(Y),1)
  feature_num <- length(X[1,])
  n <- length(X[,1])
  cut_feature = 3*n/log(n)
  dcov_store1 <- rep(0,feature_num)
  X <- scale(X)
  d_Y <- dist(Y)
  for (i in 1:feature_num) {
    dcov_store1[i] <- dcov(dist(X[,i]),d_Y)
  }
  sort_d <- sort(dcov_store1,decreasing = TRUE,index.return=TRUE)
  store_a <- NULL
  store_b <- NULL
  store_mu <- NULL
  store_lambda <- NULL
  store_a[1] <- 2
  store_b[1] <- min(20,floor((length(X[,1])/log(length(X[,1])))^0.5))
  store_mu[1] <- floor(store_a[1] + 0.618 * (store_b[1] - store_a[1]))
  store_lambda[1] <- floor(store_a[1] + 0.382 * (store_b[1] - store_a[1]))
  i = 1
  store_opt <- NULL
  store_opt[[min(20,floor((length(X[,1])/log(length(X[,1])))^0.5))+1]] <- 1
  value_left <- scalar_f_opt(X,Y,i=store_lambda[i],sort_d = sort_d,cut_feature = cut_feature)
  value_right <- scalar_f_opt(X,Y,i=store_mu[i],sort_d = sort_d,cut_feature = cut_feature)
  store_opt[[store_lambda[i]]] <- value_left
  store_opt[[store_mu[i]]] <- value_right
  while (TRUE) {
    if ((store_b[i] - store_a[i]) <= 1){
      break
    }
    if(i > 1){
      if(store_a[i] == store_a[i-1]){
        value_right <- value_left
        if(is.null(store_opt[[store_lambda[i]]])){
          value_left <- scalar_f_opt(X,Y,i=store_lambda[i],sort_d = sort_d,cut_feature = cut_feature)
          store_opt[[store_lambda[i]]] <- value_left
        }else{
          value_left <- store_opt[[store_lambda[i]]]
        }
      }
      if(store_b[i] == store_b[i-1]){
        value_left <- value_right
        if(is.null(store_opt[[store_mu[i]]])){
          value_right <- scalar_f_opt(X,Y,store_mu[i],sort_d = sort_d,cut_feature = cut_feature)
          store_opt[[store_mu[i]]] <- value_right
        }else{
          value_right <- store_opt[[store_lambda[i]]]
        }
      }
    }
    if(2 * store_lambda[i] - length(X[,1]) * log(abs(value_left$dcov)/length(X[,1])) > 2 * store_mu[i] - length(X[,1]) * log(abs(value_right$dcov)/length(X[,1]))){
      store_a[i+1] <- store_lambda[i]
      store_b[i+1] <- store_b[i]
      store_lambda[i+1] <- store_mu[i]
      store_mu[i+1] <- floor(store_a[i+1] + 0.618*(store_b[i+1] - store_a[i+1]))
    }else{
      store_a[i+1] <- store_a[i]
      store_b[i+1] <- store_mu[i]
      store_mu[i+1] <- store_lambda[i]
      store_lambda[i+1] <- floor(store_a[i+1] + 0.382*(store_b[i+1] - store_a[i+1]))
    }
    print(c(store_a[i+1],store_b[i+1]))
    i <- i +1
  }
  store_aic <- NULL
  for (i in 1:(length(store_opt)-1)) {
    if (is.null(store_opt[[i]])) store_aic[i] <- NA
    else store_aic[i] <- 2 * i - length(X[,1]) * log(abs(store_opt[[i]]$dcov)/length(X[,1]))
  }
  return(store_opt[[which.min(store_aic)]])
}


#### scalar_ifpls_test_part
#### the whole test part
# gamma_test <- function(store_dcov,store_dcov_ini,n){
#   mean_boot <- mean(n*store_dcov)
#   var_boot <- var(n*store_dcov)
#   gamma_scale <- var_boot / mean_boot
#   gamma_shape <- mean_boot / gamma_scale
#   return(1-pgamma(n*store_dcov_ini,shape = gamma_shape,  scale = gamma_scale))
# }

scalar_ifpls_test <- function(X,img, per.num = 50){
  X <- scale(X)
  dcov.ini <- scalar_fpls_select(X,Y)
  dcov.per <- NULL
  n <- dim(X)[1]
  cut_feature = 3*n/log(n)
  sel.num <- length(dcov.ini$index)
  for (i in seq(per.num)) {
    sample_index <- sample(n,n)
    dcov_store1 <- rep(0,feature_num)
    for (j in 1:feature_num) {
      dcov_store1[j] <- dcovU(X[sample_index,j],Y)
    }
    sort_d <- sort(dcov_store1,decreasing = TRUE,index.return=TRUE)
    d.per <- scalar_f_opt(X[sample_index,],Y,sel.num,sort_d = sort_d,cut_feature = cut_feature)
    dcov.per[i] <- dcov(X[sample_index,d.per$index] %*% d.per$X_par,Y)
  }
  p.value <- gamma_test(dcov.per^2,dcov(X[,dcov.ini$index] %*% dcov.ini$X_par,Y)^2,n)
  return(list(dcovU = dcov.ini$dcov, p.value = p.value))
}

#####
# #### ifpls test part
# local_opt_single <- function(X,Y,index){
#   n <- length(as.matrix(X)[,1])
#   feature_num <- length(as.matrix(X)[1,])
#   X <- scale(X)
#   ###### Determines the initial parameter values
#   sqp_par <- matrix(runif(select_num-1,0,pi))
#   cat("initial par:",beta_vector(sqp_par),"\n")
#   cat("initial dcov:",dcovU(X[,index]%*% beta_vector(sqp_par), Y),"\n")
#   if (length(X[,1]) > 2000){
#     sample_n <- sample(length(X[,1]),2000)
#     X_sub <- X[sample_n,]
#     Y_sub <- Y[sample_n,]
#   }else{
#     X_sub <- X
#     Y_sub <- Y
#   }
#   new_index <- index
#   new_beta <- sqp_par
#
#   run1_opt <- scalar_run1(40,X,Y,X_sub, Y_sub,index,sqp_par)
#   # sqp_par_new <- run1_opt$min_beta
#   cat("optimize dcov:",dcovU(X[,index]%*% beta_vector(sqp_par_new), Y),"\n")
#   return(list(dcovU = dcovU(X[,index]%*%beta_vector(sqp_par), img %*% b_hat),X_par = beta_vector(run1_opt$min_beta)))
# }
# local_permute_single <- function(time,f_index,X,Y) {
#   n <- dim(X)[1]
#   sample_index  <- sample(n,n)
#   f_mar <- local_opt_single(X[sample_index,],Y,f_index)
#   return(list(sample_index = sample_index, f_mar = f_mar))
# }
# single_ifpls_test <- function(X,Y,index.num = -1, per.num = 50,fpls_index = -1){
#   n <- dim(X)[1]
#   if(fpls_index == -1){
#     ifpls_opt <- scalar_ifpls_select(X,Y)
#     fpls_index <- ifpls_opt$index
#   }else{
#     ifpls_opt <- local_opt_single(X,Y,fpls_index)
#   }
#   permute_opt <- list()
#   for (time in 1:per.num) {
#     permute_opt[[time]] <- local_permute_single(time,fpls_index,X,Y)
#   }
#   if(index.num == -1){
#     p_store <- rep(-1,dim(X)[2])
#     for (i in 1:dim(X)[2]) {
#       if(i %in% fpls_index){
#         dcov_per <- NULL
#         for (j in 1:per.num) {
#           dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,i],Y)
#         }
#         p_store[i] <- gamma_test(dcov_per^2,dcov(X[,i],Y)^2,n)
#       }else{
#         dcov_per <- NULL
#         for (j in 1:per.num) {
#           dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,i],Y)
#         }
#         p_store[i] <- gamma_test(dcov_per^2,dcov(X[,i],Y)^2,n)
#       }
#     }
#   }else{
#     if(index.num %in% fpls_index){
#       dcov_per <- NULL
#       for (j in 1:per.num) {
#         dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,index.num],img%*%permute_opt[[j]]$f_mar$bhat)
#       }
#       p_store <- gamma_test(dcov_per^2,dcov(X[,index.num],img%*%ifpls_opt$bhat)^2,n)
#     }else{
#       dcov_per <- NULL
#       for (j in 1:per.num) {
#         dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,index.num],img%*%ifpls_opt$bhat)
#       }
#       p_store <- gamma_test(dcov_per^2,dcov(X[,index.num],img%*%ifpls_opt$bhat)^2,n)
#     }
#   }
#   return(list(dcovU = -dcov.ini$img_par$min_obj, p.value = p_store))
# }
#######

# generate_data <- function(kk){
#   n <- 200
#   p <- 500
#   feature_num = p
#   sigma_sqrt <- matrix(NA,feature_num,feature_num)
#   for (i in seq(feature_num)) {
#     for (j in seq(feature_num)) {
#       if (i == j){
#         sigma_sqrt[i,j] = 1
#       }
#       else
#         sigma_sqrt[i,j] = 0.5^abs(i-j)
#     }
#   }
#   E <- eigen(sigma_sqrt)
#   V <- E$vectors
#   sigma_matrix <- V %*% diag(sqrt(E$values)) %*% t(V)
#   X <- matrix(rnorm(n*p),n,p)
#   X <- X%*%sigma_matrix
#
#
#   img_size <- 200
#   sigma_sqrt_eta <- matrix(NA,img_size,img_size)
#   for (i in seq(img_size)) {
#     for (j in seq(img_size)) {
#       if (i == j){
#         sigma_sqrt_eta[i,j] = 1
#       }
#       else
#         sigma_sqrt_eta[i,j] = 0.5^abs(i-j)
#     }
#   }
#   E <- eigen(sigma_sqrt_eta)
#   V <- E$vectors
#   sigma_matrix <- V %*% diag(sqrt(E$values)) %*% t(V)
#
#   eta_s <- matrix(rnorm(n*img_size),n,img_size)
#   eta_s <- eta_s%*%sigma_matrix
#
#   epslion <- matrix(rnorm(n*img_size),n,img_size)
#   img <- matrix(0,n,img_size)
#   for (i in 1:n) {
#     img[i,] <- eta_s[i,]  + epslion[i,]
#     for (j in 1:img_size) {
#       if(j>50 & j <= 150){
#         img[i,j] <-  img[i,j] + ((2/sqrt(13))*((0.002*j)^(-1))*X[i,1] + (2/sqrt(13))*((0.002*j)^(-1))*(X[i,2])+
#                                    (2/sqrt(13))*((0.002*j)^(-1))*X[i,3] + (1/sqrt(13))*((0.002*j)^(-1))*X[i,30])^2
#       }
#     }
#   }
#
#   return(list(snps = X, img = img))
# }
# gene_data <- generate_data(1)
# t0 <- Sys.time()
# X <- gene_data$snps
# Y <- gene_data$img[,100]
# sc_f <- scalar_fpls_select(X,Y)
# Sys.time() - t0
