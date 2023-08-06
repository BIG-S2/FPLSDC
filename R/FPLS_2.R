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
obj_sqp <- function(beta,index,X,img,b_hat){
  return(-dcovU(X[,index] %*% beta, img%*%b_hat))
}
cons <- function(x){
  return(t(x) %*% x - 1)
}
metrospolis <- function(f,f_new,T){
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
grad1 <- function(beta,X_sub,img_sub,K_b,gamma0,index,n){
  -grad_dcov_U(sign(B_distance_2(X_sub[,index] %*% beta_vector(beta))), B_distance(img_sub %*% K_b %*% beta_vector(gamma0)), X_sub[,index]%*% grad_beta_vector(beta),rowSumsC(B_distance(img_sub %*% K_b %*% beta_vector(gamma0))),sum(rowSumsC(B_distance(img_sub %*% K_b %*% beta_vector(gamma0)))),min(n,2000))
}
generate1 <- function(beta1,X_sub,img_sub,K_b,gamma0,index,n){
  beta_new <- beta1 - runif(1,0,pi) * grad1(beta1,X_sub,img_sub,K_b,gamma0,index,n)
  return(beta_new)
}
generate1_rand <- function(beta1,beta2,X_sub,img_sub,K_b,gamma0,index,n){
  rand <- runif(length(beta1),0,pi)
  beta_new <- beta1 - runif(1,0,pi) * (grad1(beta1,X_sub,img_sub,K_b,gamma0,index,n) + sign(grad1(beta2,X_sub,img_sub,K_b,gamma0,index,n))*rand)
  return(beta_new)
}
run1 <- function(times,beta,index,X,img,X_sub,img_sub,b_hat,gamma0,K_b,n){
  f_list <- NULL
  beta_count <- NULL
  counts <- 2
  f_list[1] <- obj_sqp(beta_vector(beta),index,X,img,b_hat)
  beta_count[[1]] <- beta
  breaknum <- 0
  beta_new <- generate1(beta,X_sub,img_sub,K_b,gamma0,index,n)
  beta1 <- beta
  T <- 10
  for (i in 1:times) {
    f <- obj_sqp(beta_vector(beta),index,X,img,b_hat)
    f_new <- obj_sqp(beta_vector(beta_new),index,X,img,b_hat)
    m_num <- metrospolis(f = f,f_new = f_new,T)
    if(m_num == 0){
      breaknum <- breaknum + 1
      if(breaknum > 5) {
        break
      }
      beta_new <- generate1_rand(beta,beta1,X_sub,img_sub,K_b,gamma0,index,n)
    }
    if(m_num == 1){
      beta_count[[counts]] <- beta_new
      f_list[counts] <- f_new
      beta1 <- beta
      beta <- beta_new
      counts <- counts + 1
      breaknum <- 1
      if(f > f_new){
        beta_new <- generate1(beta,X_sub,img_sub,K_b,gamma0,index,n)
      }else{
        beta_new <- generate1_rand(beta,beta1,X_sub,img_sub,K_b,gamma0,index,n)
      }
    }
    T <- T * 0.9
  }
  min_index <- which.min(f_list)
  return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
}
#######################
# the second opt
grad2 <- function(beta,index,X_sub,img_sub,K_b,sqp_par_new,n){
  -grad_dcov_U(sign(B_distance_2(img_sub %*% K_b %*% beta_vector(beta))), B_distance(X_sub[,index] %*%beta_vector(sqp_par_new) ), img_sub %*% K_b %*% grad_beta_vector(beta),rowSumsC(B_distance(X_sub[,index] %*% beta_vector(sqp_par_new))),sum(rowSumsC(B_distance(X_sub[,index] %*% beta_vector(sqp_par_new)))),min(n,2000))
}
obj_sqp_basis <- function(sqp_beta0,X,Z,new_index,new_beta,n){
  return(-dcovU(X[,new_index] %*% beta_vector(new_beta), Z%*%sqp_beta0))
}
generate2 <- function(beta1,index,X_sub,img_sub,K_b,sqp_par_new,n){
  #step <- runif(1,0,10)
  beta_new <- beta1 - runif(1,0,pi) * grad2(beta1,index,X_sub,img_sub,K_b,sqp_par_new,n)
  return(beta_new)
}
generate2_rand <- function(beta1,beta2,index,X_sub,img_sub,K_b,sqp_par_new,n){
  rand <- runif(length(beta1),0,pi)
  beta_new <- beta1 - runif(1,0,pi) * (grad2(beta1,index,X_sub,img_sub,K_b,sqp_par_new,n) + sign(grad2(beta2,index,X_sub,img_sub,K_b,sqp_par_new,n))*rand)
  return(beta_new)
}
run2 <- function(times,beta,index,X,Z,X_sub,img_sub,new_index,new_beta,sqp_par_new,K_b,n){
  f_list <- NULL
  beta_count <- NULL
  counts <- 2
  f_list[1] <- obj_sqp_basis(beta_vector(beta),X,Z,new_index,new_beta)
  beta_count[[1]] <- beta
  breaknum <- 0
  beta_new <- generate2(beta,index,X_sub,img_sub,K_b,sqp_par_new,n)
  beta1 <- beta
  T = 10
  for (i in 1:times) {
    f <- obj_sqp_basis(beta_vector(beta),X,Z,new_index,new_beta)
    f_new <- obj_sqp_basis(beta_vector(beta_new),X,Z,new_index,new_beta)
    m_num <- metrospolis(f = f,f_new = f_new,T)
    if(m_num == 0){
      breaknum <- breaknum + 1
      if(breaknum > 5) {
        break
      }
      beta_new <- generate2_rand(beta,beta1,index,X_sub,img_sub,K_b,sqp_par_new,n)
    }
    if(m_num == 1){
      beta_count[[counts]] <- beta_new
      f_list[counts] <- f_new
      beta1 <- beta
      beta <- beta_new
      counts <- counts + 1
      breaknum <- 1
      if(f > f_new){
        beta_new <- generate2(beta,index,X_sub,img_sub,K_b,sqp_par_new,n)
      }else{
        beta_new <- generate2_rand(beta,beta1,index,X_sub,img_sub,K_b,sqp_par_new,n)
      }
    }
    T <- T * 0.95
  }
  min_index <- which.min(f_list)
  return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
}

###########################

# the third opt
grad3 <- function(beta,X_sub,new_index,img_sub ,b_hat,n){
  return(-grad_dcov_U(sign(B_distance_2(X_sub[,new_index] %*% beta_vector(beta))), B_distance(img_sub %*% b_hat), X_sub[,new_index]%*% grad_beta_vector(beta),rowSumsC(B_distance(img_sub %*% b_hat)),sum(rowSumsC(B_distance(img_sub %*% b_hat))),min(n,2000)))
}
obj_sqp_newbeta <- function(beta,new_index,X,img,b_hat){
  return(-dcovU(X[,new_index] %*% beta, img%*%b_hat))
}
generate3 <- function(beta1,X_sub,new_index,img_sub ,b_hat,n){
  beta_new <- beta1 - runif(1,0,pi) * grad3(beta1,X_sub,new_index,img_sub ,b_hat,n)
  return(beta_new)
}
generate3_rand <- function(beta1,beta2,X_sub,new_index,img_sub ,b_hat,n){
  rand <- runif(length(beta1),0,pi)
  beta_new <- beta1 - runif(1,0,pi) * (grad3(beta1,X_sub,new_index,img_sub ,b_hat,n) + sign(grad3(beta2,X_sub,new_index,img_sub ,b_hat,n))*rand)
  return(beta_new)
}
run3 <- function(times,beta,new_index,X,img,X_sub,img_sub,b_hat,n){
  f_list <- NULL
  #T <- 100
  beta_count <- NULL
  counts <- 2
  f_list[1] <- obj_sqp_newbeta(beta_vector(beta),new_index,X,img,b_hat)
  beta_count[[1]] <- beta
  breaknum <- 0
  beta_new <- generate3(beta,X_sub,new_index,img_sub,b_hat,n)
  beta1 <- beta
  T = 10
  for (i in 1:times) {
    f <- obj_sqp_newbeta(beta_vector(beta),new_index,X,img,b_hat)
    f_new <- obj_sqp_newbeta(beta_vector(beta_new),new_index,X,img,b_hat)
    m_num <- metrospolis(f = f,f_new = f_new,T)
    if(m_num == 0){
      breaknum <- breaknum + 1
      if(breaknum > 5) {
        break
      }
      beta_new <- generate3_rand(beta,beta1,X_sub,new_index,img_sub,b_hat,n)
    }
    if(m_num == 1){
      beta_count[[counts]] <- beta_new
      f_list[counts] <- f_new
      beta1 <- beta
      beta <- beta_new
      counts <- counts + 1
      breaknum <- 1
      if(f > f_new){
        beta_new <- generate3(beta,X_sub,new_index,img_sub,b_hat,n)
      }else{
        beta_new <- generate3_rand(beta,beta1,X_sub,new_index,img_sub,b_hat,n)
      }
    }
    T <- T * 0.95
  }
  min_index <- which.min(f_list)
  return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
}

# the forth opt
grad4 <- function(beta,Z_sub,X_sub,new_index,new_beta,n){
  return(-grad_dcov_U(sign(B_distance_2(Z_sub%*%beta_vector(beta))), B_distance(X_sub[,new_index] %*% beta_vector(new_beta)), Z_sub%*%grad_beta_vector(beta),rowSumsC(B_distance(X_sub[,new_index] %*% beta_vector(new_beta))),sum(rowSumsC(B_distance(X_sub[,new_index] %*% beta_vector(new_beta)))),min(n,2000)))
}
generate4 <- function(beta1,Z_sub,X_sub,new_index,new_beta,n){
  beta_new <- beta1 - runif(1,0,pi) * grad4(beta1,Z_sub,X_sub,new_index,new_beta,n)
  return(beta_new)
}
generate4_rand <- function(beta1,beta2,Z_sub,X_sub,new_index,new_beta,n){
  rand <- runif(length(beta1),0,pi)
  beta_new <- beta1 - runif(1,0,pi) * (grad4(beta1,Z_sub,X_sub,new_index,new_beta,n) + sign(grad4(beta2,Z_sub,X_sub,new_index,new_beta,n))*rand)
  return(beta_new)
}
obj_sqp_basis <- function(sqp_beta0,X,Z,new_index,new_beta,n){
  return(-dcovU(X[,new_index] %*% beta_vector(new_beta), Z%*%sqp_beta0))
}
run4 <- function(times,beta,X,Z,Z_sub,X_sub,new_index,new_beta,n){
  f_list <- NULL
  beta_count <- NULL
  counts <- 2
  f_list[1] <- obj_sqp_basis(beta_vector(beta),X,Z,new_index,new_beta)
  beta_count[[1]] <- beta
  breaknum <- 0
  beta_new <- generate4(beta,Z_sub,X_sub,new_index,new_beta,n)
  beta1 <- beta
  T = 10
  for (i in 1:times) {
    f <- obj_sqp_basis(beta_vector(beta),X,Z,new_index,new_beta)
    f_new <- obj_sqp_basis(beta_vector(beta_new),X,Z,new_index,new_beta)
    m_num <- metrospolis(f = f,f_new = f_new,T)
    if(m_num == 0){
      breaknum <- breaknum + 1
      if(breaknum > 5) {
        break
      }
      beta_new <- generate4_rand(beta,beta1,Z_sub,X_sub,new_index,new_beta,n)
    }
    if(m_num == 1){
      beta_count[[counts]] <- beta_new
      f_list[counts] <- f_new
      beta1 <- beta
      beta <- beta_new
      counts <- counts + 1
      breaknum <- 1
      if(f > f_new){
        beta_new <- generate4(beta,Z_sub,X_sub,new_index,new_beta,n)
      }else{
        beta_new <- generate4_rand(beta,beta1,Z_sub,X_sub,new_index,new_beta,n)
      }
    }
    T <- T *0.95
  }
  min_index <- which.min(f_list)
  return(list(break_num = breaknum,T = T, min_beta = beta_count[[min_index]],min_value = f_list[min_index], beta = beta_count, f_list = f_list))
}
######## function
norm_1 <- function(pro_par){
  return(pro_par/sqrt(as.numeric(t(pro_par)%*% pro_par)))
}
obj_BB <-function(theta,X,index,pro_par,sort_index,select_num,add_i,img, b_hat){
  return(-dcovU(cos(theta)*(X[,index] %*% norm_1(pro_par))+sin(theta)*X[,sort_index$ix[-c(1:select_num)][add_i]],img %*% b_hat))
}
obj_grad_BB <- function(theta,X,index,pro_par,sort_index,select_num,add_i,img,b_hat){
  n <- length(X[,1])
  if (length(X[,1]) > 2000){
    sample_index <- sample(length(X[,1]),2000)
    grad_ratio <- -grad_dcovU_ratio(sign(B_distance_2(cos(theta)*(X[sample_index,index] %*% norm_1(pro_par))+sin(theta)*X[sample_index,sort_index$ix[-c(1:select_num)][add_i]])),
                                    B_distance_2(-sin(theta)*(X[sample_index,index] %*% norm_1(pro_par))+cos(theta)*X[sample_index,sort_index$ix[-c(1:select_num)][add_i]]),
                                    B_distance((img %*% b_hat)[sample_index,]),rowSumsC(B_distance((img %*% b_hat)[sample_index,])),sum(rowSumsC(B_distance((img %*% b_hat)[sample_index,]))),length(sample_index))
  }
  else{
    grad_ratio <- -grad_dcovU_ratio(sign(B_distance_2(cos(theta)*(X[,index] %*% norm_1(pro_par))+sin(theta)*X[,sort_index$ix[-c(1:select_num)][add_i]])),
                                    B_distance_2(-sin(theta)*(X[,index] %*% norm_1(pro_par))+cos(theta)*X[,sort_index$ix[-c(1:select_num)][add_i]]),
                                    B_distance(img %*% b_hat),rowSumsC(B_distance(img %*% b_hat)),sum(rowSumsC(B_distance(img %*% b_hat))),n)
  }
  return(grad_ratio)
}
grad_norm <- function(x){
  return(sqrt(t(x) %*% x))
}
BB_opt <- function(x0,X,index,pro_par,sort_index,select_num,add_i,img, b_hat){
  x0 <- 0.01
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
  store_value[1] <- obj_BB(store_x[1],X,index,pro_par,sort_index,select_num,add_i,img, b_hat)
  #t0 = Sys.time()
  while (TRUE) {
    for (i in 1:100) {
      if(obj_BB(store_x[k] - alpha*obj_grad_BB(store_x[k],X,index,pro_par,sort_index,select_num,add_i,img,b_hat),X,index,pro_par,sort_index,select_num,add_i,img, b_hat) < max(store_value[(k + 1 - min(k,M)):k]) - c1*alpha*grad_norm(obj_grad_BB(store_x[k],X,index,pro_par,sort_index,select_num,add_i,img,b_hat))^2){
        break
      }
      alpha <- beta * alpha
    }
    store_x[k+1] <- store_x[k] - alpha*t(obj_grad_BB(store_x[k],X,index,pro_par,sort_index,select_num,add_i,img,b_hat))
    s <- store_x[k+1] - store_x[k]
    y <- obj_grad_BB(store_x[k+1],X,index,pro_par,sort_index,select_num,add_i,img,b_hat) - obj_grad_BB(store_x[k],X,index,pro_par,sort_index,select_num,add_i,img,b_hat)
    alpha <- as.numeric((t(s) %*% y)/(t(y) %*% y))
    if(is.null(alpha)) alpha <- 0.01
    if(alpha < 0.01) alpha <- 0.01
    if(alpha > 0.5) alpha <- 0.5
    k <- k+1
    if(grad_norm(alpha * obj_grad_BB(store_x[k],X,index,pro_par,sort_index,select_num,add_i,img,b_hat)) < epslion) break
    store_value[k] <- obj_BB(store_x[k],X,index,pro_par,sort_index,select_num,add_i,img, b_hat)
    if(abs(store_value[k] - store_value[(k-1)]) < 1e-3) break
    if(k >= max_k) break
  }
  return(list(par = store_x, value = -store_value))
}

# APLS select the basis function
APLS <- function(x, Y, p){
  n = length(x)
  d = dim(Y)[2]
  K_b = matrix(0, d, p)
  K_b[,1] = (t(Y) - apply(Y, 2, mean))%*%(x - mean(x))/n
  if(p > 1){
    for(i in 2:p){
      K_b[,i] = (t(Y) - apply(Y, 2, mean))%*%(t(t(Y) - apply(Y, 2, mean))%*%K_b[,i-1])/n
    }
  }
  return(K_b)
}
# orthogonalization of basis function
orth <- function(basis){
  n = dim(basis)[1]
  p = dim(basis)[2]
  ortho_basis = matrix(0,n,p)
  ortho_basis[,1] = basis[,1]/sqrt(sum(basis[,1]^2))
  if(p>1){
    for(i in 2:p){
      ortho_basis[,i] = basis[,i]
      for(j in 1:(i-1)){
        ortho_basis[,i] = ortho_basis[,i] - sum(basis[,i]*ortho_basis[,j])*ortho_basis[,j]
      }
      ortho_basis[,i] = ortho_basis[,i]/sqrt(sum(ortho_basis[,i]^2))
    }
  }
  return(ortho_basis)
}
v_norm <- function(x){
  return(x/(as.vector(t(x) %*% x)**0.5))
}

best_select <- function(X,img,sort_d,select_num = 4,p = 5,cut_feature,true_index = -1,min_feature = 5){
  n <- length(as.matrix(X)[,1])
  feature_num <- length(as.matrix(X)[1,])
  sort_index <- sort_d
  # feature_num_sub <- floor(min(5*n/log(n),length(as.matrix(X)[1,])))
  feature_num_sub <- floor(min(cut_feature,length(as.matrix(X)[1,])))
  sort_index_sub <- list(x = sort_d$x[1:feature_num_sub],ix = sort_d$ix[1:feature_num_sub])
  index <- sort_index_sub$ix[1:select_num]
  cat("the initial index:",index,"\n")
  # cat("\n")
  if (true_index[1] != -1) cat("Before the update, the position of the active variable：",match(true_index,sort_d$ix),"\n")

  ###### Determines the initial parameter values
  sqp_par <- matrix(runif(select_num-1,0,pi))
  gamma0 <- matrix(runif(p-1,0,pi))
  K_b <- APLS(X[,index] %*%beta_vector(sqp_par), img, p)
  K_b <- orth(K_b)
  b_hat <- K_b %*% beta_vector(gamma0)
  # cat("initial dcov:",dcovU(X[,index]%*% beta_vector(sqp_par), img %*% b_hat),"\n")
  # Iteratively optimize initial parameters
  if (length(X[,1]) > 2000){
    sample_n <- sample(length(X[,1]),2000)
    X_sub <- X[sample_n,]
    img_sub <- img[sample_n,]
  }else{
    X_sub <- X
    img_sub <- img
  }
  new_index <- index
  new_beta <- sqp_par
  ini.times <- 1
  while (TRUE) {
    dcov_old <- dcovU(X[,index]%*%beta_vector(sqp_par), img %*% K_b %*%beta_vector(gamma0))
    sqp_par_new <- run1(40,sqp_par,index,X,img,X_sub,img_sub,b_hat,gamma0,K_b,n)$min_beta
    K_b <- APLS(X[,index] %*% beta_vector(sqp_par_new), img, p)
    K_b <- orth(K_b)
    Z <- img %*% K_b
    gamma0_new <- run2(40,gamma0,index,X,Z,X_sub,img_sub,new_index,new_beta,sqp_par_new,K_b,n)$min_beta
    dcov_new <- dcovU(X[,index]%*% beta_vector(sqp_par_new), img %*% K_b %*% beta_vector(gamma0_new))
    # cat(dcov_new,dcov_old,"\n")
    if (dcov_new <= dcov_old){
      break
    }else{
      gamma0 <- gamma0_new
      sqp_par <- sqp_par_new
    }
    ini.times <- ini.times + 1
    if(ini.times > 10){break}
  }

  # cat("optimate initial dcov:",dcov_old,"\n")
  K_b <- APLS(X[,index] %*% beta_vector(sqp_par), img, p)
  K_b <- orth(K_b)
  b_hat <- K_b %*% beta_vector(gamma0)
  # dcovU(X[,index]%*% beta_vector(sqp_par_new), img %*% b_hat)
  ###### store the variation of dcovU of each variable
  iterations <- 1
  end_dcov <- NULL  # store the iterated dcov
  end_dcov[1] <- 0
  while (TRUE) {
    index_variation <- rep(-1,feature_num)
    d_cov <- dcovU(X[,index]%*%beta_vector(sqp_par),(img %*% b_hat))
    sqp_par1 <- beta_vector(sqp_par)
    for (min_i in 1:length(index)) {
      pro_par <- sqp_par1[-min_i]/sqrt(sum((sqp_par1[-min_i])^2))
      index_variation[index[min_i]] <-  abs(d_cov - dcovU(as.matrix(X[,index[-min_i]])%*%pro_par,(img %*% b_hat)))
    }
    pro_par <- sqp_par1 / sqp_par1[1]
    add_length <- max((select_num+1),feature_num_sub)
    for (add_i in 1:length(seq(add_length-select_num))) {
      # cat(add_i)
      index_variation[sort_index$ix[-c(1:select_num)][add_i]] <- abs(tail(BB_opt(0.01,X,index,pro_par,sort_index,select_num,add_i,img, b_hat)$value,1) - d_cov)
    }
    sort_index_sub <- sort(index_variation,decreasing = TRUE,index.return=TRUE)
    # sort_index_sub <- list(x = c(sort_index_sub2$x[1:feature_num_sub],sort_d$x[feature_num_sub]), ix = sort_index_sub2$ix[1:feature_num_sub])
    sort_index_old <- sort_index
    sort_index$ix[1:feature_num_sub] <- sort_index_sub$ix[1:feature_num_sub]

    new_index <- sort_index$ix[1:select_num]

    cat("the new updated index:",new_index,"\n")
    if (true_index[1] != -1) cat("After the update, the position of the active variable：", match(true_index,sort_index$ix),"\n")

    ##### After the new index is determined, Update the coefficient of new vaiable
    if (n > 2000){
      sample_n <- sample(n,2000)
      X_sub <- X[sample_n,]
      img_sub <- img[sample_n,]
    }else {
      X_sub <- X
      img_sub <- img
    }
    ini_beta <- runif(length(new_index)-1,0,pi)
    new_beta <- run3(40,ini_beta,new_index,X,img,X_sub,img_sub,b_hat,n)$min_beta
    # cat("update the coefficient of  new index, the new dcovU ",dcovU(X[,new_index] %*%beta_vector(new_beta) , img %*% b_hat),"\n")

    b_hat_old <- b_hat
    K_b <- APLS((X[,new_index] %*% beta_vector(new_beta)) , img, p)
    K_b <- orth(K_b)
    Z <- (img %*% K_b)
    if (n > 2000){
      sample_n <- sample(n,2000)
      X_sub <- X[sample_n,]
      Z_sub <- Z[sample_n,]
    }else {
      X_sub <- X
      Z_sub <- Z
    }
    beta0 <- runif(p-1,0,pi)
    b_hat <- K_b %*%beta_vector(run4(40,beta0,X,Z,Z_sub,X_sub,new_index,new_beta,n)$min_beta)
    if (dcovU(X[,new_index] %*%beta_vector(new_beta) ,img %*% b_hat) <= dcovU(X[,new_index] %*%beta_vector(new_beta),img %*% b_hat_old)){
      b_hat <- b_hat_old
    }
    updata_dcov <- dcovU(X[,new_index] %*%beta_vector(new_beta) , img %*% b_hat)
    # cat("update the coefficient of basis, the upated dcovU ",updata_dcov,"\n")
    # cat("\n")
    end_dcov[iterations+1] <- updata_dcov

    if (iterations+1 > 10){
      cat("Reaches the maximum number of iterations\n")
      if(true_index[1] != -1)return(list(index= sort_index$ix[1:select_num], position = sort_index$ix, value = updata_dcov, par_x = beta_vector(new_beta), par_img = b_hat))
      else return(list(index= sort_index$ix[1:select_num],position = sort_index$ix,value = updata_dcov, par_x = beta_vector(new_beta), par_img = b_hat))
      break
    }
    if ((end_dcov[iterations]/end_dcov[iterations+1] - 1)>0.01){
      cat("The dcovU decrease, stop!\n")
      if(true_index[1] != -1)return(list(index= sort_index_old$ix[1:select_num], position = sort_index_old$ix,value = dcovU(X[,index] %*% beta_vector(sqp_par), img %*% b_hat_old),par_x =  beta_vector(sqp_par), par_img = b_hat_old))
      else return(list(index= sort_index_old$ix[1:select_num],position = sort_index_old$ix,value = dcovU(X[,index] %*%beta_vector(sqp_par), img %*% b_hat_old),par_x =  beta_vector(sqp_par), par_img = b_hat_old))
      break
    }
    if (sum(sort(index) == sort(new_index)) == select_num){
      cat("The iteration stop\n")
      if(true_index[1] != -1)return(list(index= sort_index$ix[1:select_num], position = sort_index$ix,value = updata_dcov,par_x = beta_vector(new_beta) , par_img = b_hat))
      else return(list(index= sort_index$ix[1:select_num],position = sort_index$ix,value = updata_dcov, par_x = beta_vector(new_beta), par_img = b_hat))
      break
    }
    index <- new_index
    iterations <- iterations + 1
    sqp_par <- new_beta
    feature_num_sub <- feature_num_sub/min_feature
  }
  Sys.time() - t0
}

ifplsdc_selnum <- function(X,img,sel_num,sort_d,cut_feature = NULL){
  if(is.null(cut_feature)){
    cut_feature = dim(X)[2]
  }
  p <- min(5,dim(img)[2])
  store_value <- best_select(X,img,select_num = sel_num,p= 5,cut_feature = cut_feature, min_feature = 5,sort_d = sort_d)
  cat("the selected variable:",store_value$index,"\n")
  return(list(index = store_value$index, rank = store_value$position,dcovU = store_value$value, par = store_value$par_x,bhat = store_value$par_img))
}
fpls_dc_ix <- function(X,Y){
  dcov_vec <- NULL
  X <- scale(X)
  n <- dim(X)[1]
  dim_len <- dim(X)[2]
  mar_dcov <- rep(-1,dim_len)
  if (dim_len > floor(3*n/log(n))){
    d_Y <- dist(Y)
    for (i in 1:length(X[1,])) {
      dcov_vec[i] <- dcov(dist(X[,i]),d_Y)
    }
    dcov_vec_sort <- sort(dcov_vec,decreasing = TRUE, index.return = TRUE)

    dcov_vec_sort_sub <- list(x = dcov_vec_sort$x[1:floor(3*n/log(n))],ix = dcov_vec_sort$ix[1:floor(3*n/log(n))])
    for (i in 1:floor(3*n/log(n))) {
      m_dcov <- marginal_dcov(X[,dcov_vec_sort_sub$ix[i]],Y)
      mar_dcov[dcov_vec_sort_sub$ix[i]] <- -m_dcov$fplsdc.list$min_obj
    }
    sort_d <- sort(mar_dcov,decreasing = TRUE,index.return=TRUE)

    sort_index <- list(x = c(sort_d$x[1:floor(3*n/log(n))], dcov_vec_sort$x[-seq(floor(3*n/log(n)))]),
      ix = c(sort_d$ix[1:floor(3*n/log(n))], dcov_vec_sort$ix[-seq(floor(3*n/log(n)))]))
  }else{
    for (i in 1:dim_len) {
      m_dcov <- marginal_dcov(X[,i],Y)
      mar_dcov[i] <- -m_dcov$fplsdc.list$min_obj
    }
    sort_d <- sort(mar_dcov,decreasing = TRUE,index.return=TRUE)
    sort_index <- sort_d
  }
  return(rank = sort_index)
}
fpls_select <- function(X,Y){
  # library(mFPLS)
  library(energy)
  feature_num <- length(X[1,])
  n <- length(X[,1])
  cut_feature = 3*n/log(n)
  dcov_store1 <- rep(0,feature_num)
  X <- scale(X)

  if (feature_num > floor(n/log(n))){
    d_Y <- dist(Y)
    for (i in 1:feature_num) {
      dcov_store1[i] <- dcov(dist(X[,i]),d_Y)
    }
    sort_d <- sort(dcov_store1,decreasing = TRUE,index.return=TRUE)
  }else{
    sort_d <- fpls_dc_ix(X,Y)
  }

  store_a <- NULL
  store_b <- NULL
  store_mu <- NULL
  store_lambda <- NULL
  store_a[1] <- 2
  store_b[1] <- min(20,floor((length(X[,1])/log(length(X[,1])))^0.5))
  store_mu[1] <- floor(store_a[1] + 0.618 * (store_b[1] - store_a[1]))
  store_lambda[1] <- floor(store_a[1] + 0.382 * (store_b[1] - store_a[1]))
  i = 1
  store_opt <- list()
  store_opt[[min(20,floor((length(X[,1])/log(length(X[,1])))^0.5))+1]] <- 1
  value_left <- ifplsdc_selnum(X,Y,store_lambda[i],sort_d = sort_d,cut_feature = cut_feature)
  value_right <- ifplsdc_selnum(X,Y,store_mu[i],sort_d = sort_d,cut_feature = cut_feature)
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
          value_left <- ifplsdc_selnum(X,Y,store_lambda[i],sort_d = sort_d,cut_feature = cut_feature)
          store_opt[[store_lambda[i]]] <- value_left
        }else{
          value_left <- store_opt[[store_lambda[i]]]
        }
      }
      if(store_b[i] == store_b[i-1]){
        value_left <- value_right
        if(is.null(store_opt[[store_mu[i]]])){
          value_right <- ifplsdc_selnum(X,Y,store_mu[i],sort_d = sort_d,cut_feature = cut_feature)
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
    # cat("after updated",c(store_a[i+1],store_b[i+1]))
    i <- i +1
  }
  store_aic <- NULL
  for (i in 1:(length(store_opt)-1)) {
    if (is.null(store_opt[[i]])) store_aic[i] <- NA
    else store_aic[i] <- 2 * i - length(X[,1]) * log(abs(store_opt[[i]]$dcov)/length(X[,1]))
  }
  return(store_opt[[which.min(store_aic)]])
  # return(value_left)
}

ifpls_dc <- function(X,Y){
  library(energy)
  if(sum(c(is.null(dim(Y)), (dim(Y) == 1))) == 1){
    Y <- matrix(Y)
    return(scalar_fpls_select(X,Y))
  }else{
    return(fpls_select(X,Y))
  }
}
#### marginal fpls
img_grad2 <- function(img_sub,X_sub,K_b,beta){
  return(-grad_dcov_U(sign(B_distance_2(img_sub %*% K_b %*% beta_vector(beta))), B_distance(X_sub), img_sub %*% K_b %*% grad_beta_vector(beta),rowSumsC(B_distance(X_sub)),sum(rowSumsC(B_distance(X_sub))),min(n,2000)))
}
img_obj <- function(X,Z,sqp_beta0){
  return(-dcovU(X, Z%*%sqp_beta0))
}
img_generate2 <- function(img_sub,X_sub,K_b,beta1){
  beta_new <- beta1 - runif(1,0,pi) * img_grad2(img_sub,X_sub,K_b,beta1)
  return(beta_new)
}
img_generate2_rand <- function(img_sub,X_sub,K_b,beta1,beta2){
  rand <- runif(length(beta1),0,pi)
  beta_new <- beta1 - runif(1,0,pi) * (img_grad2(img_sub,X_sub,K_b,beta1) + sign(img_grad2(img_sub,X_sub,K_b,beta2))*rand)
  return(beta_new)
}
img_run2 <- function(times,beta,X,Z,X_sub,img_sub,K_b){
  f_list <- NULL
  beta_count <- NULL
  counts <- 2
  f_list[1] <- img_obj(X,Z,beta_vector(beta))
  beta_count[[1]] <- beta
  breaknum <- 0
  beta_new <- img_generate2(img_sub,X_sub,K_b,beta)
  beta1 <- beta
  T = 10
  for (i in 1:times) {
    f <- img_obj(X,Z,beta_vector(beta))
    #beta_new <- generate2(beta)
    f_new <- img_obj(X,Z,beta_vector(beta_new))
    m_num <- metrospolis(f = f,f_new = f_new,T)
    if(m_num == 0){
      breaknum <- breaknum + 1
      if(breaknum > 5) {
        break
      }
      beta_new <- img_generate2_rand(img_sub,X_sub,K_b,beta,beta1)
    }
    if(m_num == 1){
      beta_count[[counts]] <- beta_new
      f_list[counts] <- f_new
      beta1 <- beta
      beta <- beta_new
      counts <- counts + 1
      breaknum <- 1
      if(f > f_new){
        beta_new <- img_generate2(img_sub,X_sub,K_b,beta)
      }else{
        beta_new <- img_generate2_rand(img_sub,X_sub,K_b,beta,beta1)
      }
    }
    T <- T * 0.95
  }
  min_index <- which.min(f_list)
  return(list(break_num = breaknum,T = T,beta = beta_vector(beta_count[[min_index]]),min_obj= min(f_list),f_list = f_list))
}
marginal_dcov <- function(X,Y,p=5,max.time = 50){
  library(energy)
  X <- as.matrix(X)
  n <<- length(Y[,1])
  p <- min(5,dim(Y)[2])
  gamma0 <<- matrix(runif(p-1,0,pi))
  K_b <- APLS(X, Y, p)
  K_b <- orth(K_b)
  b_hat <- K_b %*% beta_vector(gamma0)
  #t0 <- Sys.time()
  # cat("initial dcov:",dcovU(X, Y %*% b_hat),"\n")
  if (length(X[,1]) > 2000){
    sample_n <- sample(length(Y[,1]),2000)
    X_sub <- X[sample_n,]
    img_sub <- Y[sample_n,]
  }else{
    X_sub <- X
    img_sub <- Y
  }
  Z <- Y %*% K_b
  gamma0_new <- img_run2(max.time,gamma0,X,Z,X_sub,img_sub,K_b)
  # cat("optimate initial dcov:",dcovU(X, img %*% K_b %*% gamma0_new$beta),"\n")
  return(list(fplsdc.list = gamma0_new, coef.fun = K_b %*% gamma0_new$beta))
}
### return the order of the variable
fpls_dc <- function(X,Y){
  dcov_vec <- NULL
  X <- scale(X)
  n <- dim(X)[1]
  dim_len <- dim(X)[2]
  mar_dcov <- rep(-1,dim_len)
  if (dim_len > floor(3*n/log(n))){
    d_Y <- dist(Y)
    for (i in 1:length(X[1,])) {
      dcov_vec[i] <- dcov(dist(X[,i]),d_Y)
    }
    dcov_vec_sort <- sort(dcov_vec,decreasing = TRUE, index.return = TRUE)

    dcov_vec_sort_sub <- list(x = dcov_vec_sort$x[1:floor(3*n/log(n))],ix = dcov_vec_sort$ix[1:floor(3*n/log(n))])
    for (i in 1:floor(3*n/log(n))) {
      m_dcov <- marginal_dcov(X[,dcov_vec_sort_sub$ix[i]],Y)
      mar_dcov[dcov_vec_sort_sub$ix[i]] <- -m_dcov$fplsdc.list$min_obj
    }
    sort_d <- sort(mar_dcov,decreasing = TRUE,index.return=TRUE)

    sort_index <- c(sort_d$ix[1:floor(3*n/log(n))], dcov_vec_sort$ix[-seq(floor(3*n/log(n)))])
  }else{
    for (i in 1:dim_len) {
      m_dcov <- marginal_dcov(X[,i],Y)
      mar_dcov[i] <- -m_dcov$fplsdc.list$min_obj
    }
    sort_d <- sort(mar_dcov,decreasing = TRUE,index.return=TRUE)
    sort_index <- sort_d$ix
  }
  return(rank = sort_index)
}
#### mfpls test part
gamma_test <- function(store_dcov,store_dcov_ini,n){
  mean_boot <- mean(n*store_dcov^2)
  var_boot <- var(n*store_dcov^2)
  gamma_scale <- var_boot / mean_boot
  gamma_shape <- mean_boot / gamma_scale
  return(1-pgamma(n*store_dcov_ini^2,shape = gamma_shape,  scale = gamma_scale))
}

fplsdc_test <- function(X,Y,index.num = NULL,per.num = 50){
  n <- dim(X)[1]
  dim_len <- dim(X)[2]
  if (is.null(index.num)){
    p.value  <- NULL
    for (k in 1:dim_len) {
      m_dcov <- marginal_dcov(X[,k],Y)
      # mar_dcov[i] <- -m_dcov$fplsdc.list$min_obj
      for (i in 1:per.num) {
        sample_index <- sample(n,n)
        dcov.per <- marginal_dcov(X[sample_index,k],Y)
        store_dcov.per[i] <- dcov(X[sample_index,k],Y%*%dcov.per$coef.fun)
      }
      p.value[k] = gamma_test(store_dcov.per,dcov(X[,index.num],Y%*%dcov.ini$coef.fun),n)
    }
    return(p.value = p.value)
  }else{
    p.value  <- NULL
    for (k in seq(length(index.num))) {
      dcov.ini <- marginal_dcov(X[,index.num[k]],Y)
      store_dcov.per <- c()
      ## normalize the X (SNP1S data), the column of X represent SNP.
      for (i in 1:per.num) {
        sample_index <- sample(n,n)
        dcov.per <- marginal_dcov(X[sample_index,index.num[k]],Y)
        store_dcov.per[i] <- dcov(X[sample_index,index.num[k]],Y%*%dcov.per$coef.fun)
      }
      p.value[k] = gamma_test(store_dcov.per,dcov(X[,index.num[k]],Y%*%dcov.ini$coef.fun),n)
    }
    return(p.value = p.value)
    }

}

gamma_test_ratio <- function(mean_boot,var_boot,store_dcov_ini,n){
  gamma_scale <- var_boot / mean_boot
  gamma_shape <- mean_boot / gamma_scale
  return(1-pgamma(n*store_dcov_ini^2,shape = gamma_shape,  scale = gamma_scale))
}

fplsdc_test_gwas <- function(X,Y,index.num = NULL,per.num = 50){
    n <- dim(X)[1]
    dim_x <- dim(X)[2]
    X1 <- as.matrix(sample(0:2, size = n, replace = TRUE, prob = c(1/3, 1/3, 1/3)))
    X1 <- scale(X1)
    store_dcov.per <- c()
    for (i in 1:per.num) {
      sample_index <- sample(n,n)
      dcov.per <- marginal_dcov(X1[sample_index],Y)
      store_dcov.per[i] <- dcov(X1[sample_index],Y%*%dcov.per$coef.fun)
    }
    mean_x <- mean(n*store_dcov.per^2)
    var_x <- var(n*store_dcov.per^2)

    dcovu.ini <- marginal_dcov(X1,Y)
    mean_dist_x <- sum(rowSums(as.matrix(dist(X1))))
    var_x_base <- dcov(X1,X1)^2

    p.value <- c()
    X <- scale(X)

    if(is.null(index.num) == TRUE){
      for (k in seq(dim_x)) {
        dcovu.ini <- marginal_dcov(X[,k],Y)
        dcov.ini <- dcov(X[,k],Y%*%dcovu.ini$coef.fun)
        mean_dist_xi <- sum(rowSums(as.matrix(dist(X[,k]))))
        mean_ratio <- ((mean_dist_xi)/(mean_dist_x))*mean_x
        var_ratio <- ((dcov(X[,k],X[,k])^2)/(var_x_base))*var_x
        p.value[k] <- gamma_test_ratio(mean_ratio,var_ratio,dcov.ini,n)
      }
    }else{
      for (k in seq(length(index.num))) {
        dcovu.ini <- marginal_dcov(X[,index.num[k]],Y)
        dcov.ini <- dcov(X[,index.num[k]],Y%*%dcovu.ini$coef.fun)
        mean_dist_xi <- sum(rowSums(as.matrix(dist(X[,index.num[k]]))))
        mean_ratio <- ((mean_dist_xi)/(mean_dist_x))*mean_x
        var_ratio <- ((dcov(X[,index.num[k]],X[,index.num[k]])^2)/(var_x_base))*var_x
        p.value[k] <- gamma_test_ratio(mean_ratio,var_ratio,dcov.ini,n)
      }
    }
  return(p.value = p.value)
}

#### ifpls test part
ifplsdc_index <- function(X,Y,index,p = 5){
  library(energy)
  feature_num <- length(X[1,])
  n <- length(X[,1])
  X <- scale(X)
  select_num <- length(index)
  p <- min(5,dim(Y)[2])
  ###### Determines the initial parameter values
  sqp_par <- matrix(runif(select_num-1,0,pi))
  gamma0 <- matrix(runif(p-1,0,pi))
  K_b <- APLS(X[,index] %*%beta_vector(sqp_par), Y, p)
  K_b <- orth(K_b)
  b_hat <- K_b %*% beta_vector(gamma0)
  # cat("initial dcov:",dcovU(X[,index]%*% beta_vector(sqp_par), img %*% b_hat),"\n")
  if (length(X[,1]) > 2000){
    sample_n <- sample(length(X[,1]),2000)
    X_sub <- X[sample_n,]
    img_sub <- Y[sample_n,]
  }else{
    X_sub <- X
    img_sub <- Y
  }
  new_index <- index
  new_beta <- sqp_par
  ini.times <- 1
  while (TRUE) {
    dcov_old <- dcovU(X[,index]%*%beta_vector(sqp_par), Y %*% K_b %*%beta_vector(gamma0))
    sqp_par_new <- run1(50,sqp_par,index,X,Y,X_sub,img_sub,b_hat,gamma0,K_b,n)$min_beta
    K_b <- APLS(X[,index] %*%beta_vector(sqp_par_new), Y, p)
    K_b <- orth(K_b)
    Z <- Y %*% K_b
    gamma0_new <- run2(50,gamma0,index,X,Z,X_sub,img_sub,new_index,new_beta,sqp_par_new,K_b,n)$min_beta
    dcov_new <- dcovU(X[,index]%*% beta_vector(sqp_par_new), Y %*% K_b %*% beta_vector(gamma0_new))
    # cat(dcov_new,dcov_old,"\n")
    if (dcov_new <= dcov_old){
      break
    }else{
      gamma0 <- gamma0_new
      sqp_par <- sqp_par_new
    }
    ini.times <- ini.times + 1
    if(ini.times > 10){break}
  }
  # cat("optimate initial dcov:",dcov_old,"\n")
  K_b <- APLS(X[,index] %*% beta_vector(sqp_par), Y, p)
  K_b <- orth(K_b)
  b_hat <- K_b %*% beta_vector(gamma0)
  return(list(dcovU = dcovU(X[,index]%*%beta_vector(sqp_par), Y %*% b_hat),par = beta_vector(sqp_par), bhat = b_hat))
}
local_permute <- function(time,f_index,X,img) {
  # library(mFPLS)
  n <- dim(X)[1]
  sample_index  <- sample(n,n)
  X1 <- X[sample_index,]
  #sample_index <- sample(n,n)
  f_mar <- ifplsdc_index(X1,img,f_index,p=5)
  # f_mar <- marginal_dcov(X1[,f$index]%*%beta_vector(f$X_par),img)
  return(list(sample_index = sample_index, f_mar = f_mar))
}

ifplsdc_test <- function(X,Y,index.num = NULL, per.num = 50,fpls_index = NULL){
  n <- dim(X)[1]
  if(is.null(fpls_index)){
    ifpls_opt <- fpls_select(X,Y)
    fpls_index <- ifpls_opt$index
    cat("\n")
    cat("The valuable variable is ",fpls_index)
    cat("\n")
  }else{
    ifpls_opt <- ifplsdc_index(X,Y,fpls_index)
  }
  cat("The permutation procedure is running...","\n")
  permute_opt <- list()
  for (time in 1:per.num) {
    permute_opt[[time]] <- local_permute(time,fpls_index,X,Y)
  }
  if(is.null(index.num)){
    p_store <- rep(-1,dim(X)[2])
    for (i in 1:dim(X)[2]) {
      if(i %in% fpls_index){
        dcov_per <- NULL
        for (j in 1:per.num) {
          dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,i],Y%*%permute_opt[[j]]$f_mar$bhat)
        }
        p_store[i] <- gamma_test(dcov_per,dcov(X[,i],Y%*%ifpls_opt$bhat),n)
      }else{
        dcov_per <- NULL
        for (j in 1:per.num) {
          dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,i],Y%*%ifpls_opt$bhat)
        }
        p_store[i] <- gamma_test(dcov_per,dcov(X[,i],Y%*%ifpls_opt$bhat),n)
      }
    }
  }else{
    p_store <- rep(-1,length(index.num))
    for (i in seq(length(index.num))) {
      if(index.num[i] %in% fpls_index){
        dcov_per <- NULL
        for (j in 1:per.num) {
          dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,index.num[i]],Y%*%permute_opt[[j]]$f_mar$bhat)
        }
        p_store[i] <- gamma_test(dcov_per,dcov(X[,index.num[i]],Y%*%ifpls_opt$bhat),n)
      }else{
        dcov_per <- NULL
        for (j in 1:per.num) {
          dcov_per[j] <- dcov(X[permute_opt[[j]]$sample_index,index.num[i]],Y%*%ifpls_opt$bhat)
        }
        p_store[i] <- gamma_test(dcov_per,dcov(X[,index.num[i]],Y%*%ifpls_opt$bhat),n)
      }
    }
  }
  return(p.value = p_store)
}

