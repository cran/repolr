`repolr` <-
function(formula,subjects,data,times,categories,
           scalevalue=1,corstr="ar1",maxiter=10,tol=0.001,
           alpha=0.5,po.test=TRUE,fixed=FALSE){

 # load functions
 
 # function to construct C matrices

 cmat.work <- function(calph,ctimes,corstr){

  # first order autoregressive model
  if(corstr=="ar1"){

  C_mat <- matrix(0,nrow=length(ctimes),ncol=length(ctimes)) 
  IC_mat <- C_mat; gIC_mat <- C_mat; ggIC_mat <- C_mat
  for (i in 1:length(ctimes)){
   for (j in 1:length(ctimes)){

   # matrix C
   diff <- abs(ctimes[i]-ctimes[j])
   C_mat[i,j] <- calph**(diff)

   # matrices IC, gIC and ggIC

   # diagonals
   if (i==1 & j==1){
    
    # top left corner
    if (calph==0){
     IC_mat[i,j] <- 1
     gIC_mat[i,j] <- 0
     ggIC_mat[i,j] <- 2
    }
    else
    {
     diff1 <- abs(ctimes[i+1]-ctimes[j])
     dd1 <- 1/(1-calph^(2*diff1))
     gdd1 <- 2*(dd1^2)*diff1*(calph^(2*diff1-1))
     ggdd1 <- 2*(dd1^3)*diff1*(calph^(2*diff1-2))*
          ((2*diff1-1)+(2*diff1+1)*calph^(2*diff1))
     IC_mat[i,j] <- dd1
     gIC_mat[i,j] <- gdd1
     ggIC_mat[i,j] <- ggdd1
    }
   }
   else if (i==length(ctimes) & j==length(ctimes))
   {
   
    # bottom right corner 
    if (calph==0){
     IC_mat[i,j] <- 1
     gIC_mat[i,j] <- 0
     ggIC_mat[i,j] <- 2
    }
    else
    {
     diff1 <- abs(ctimes[i]-ctimes[j-1])
     dd1 <- 1/(1-calph^(2*diff1))
     gdd1 <- 2*(dd1^2)*diff1*(calph^(2*diff1-1))
     ggdd1 <- 2*(dd1^3)*diff1*(calph^(2*diff1-2))*
           ((2*diff1-1)+(2*diff1+1)*calph^(2*diff1))
     IC_mat[i,j] <- dd1
     gIC_mat[i,j] <- gdd1
     ggIC_mat[i,j] <- ggdd1
    }
   }
   else if (i==j){

    # other diagonals
    if (calph==0){
     IC_mat[i,j] <- 1
     gIC_mat[i,j] <- 0
     ggIC_mat[i,j] <- 4
    }
    else
    {
     diff11 <- abs(ctimes[i]-ctimes[j-1])
     diff12 <- abs(ctimes[i+1]-ctimes[j])
     diff2 <- abs(ctimes[i+1]-ctimes[j-1])
     dd2 <- (1-calph^(2*diff2))/((1-calph^(2*diff11))*(1-calph^(2*diff12)))
     dd11 <- 1/(1-calph^(2*diff11))
     dd12 <- 1/(1-calph^(2*diff12))
     dd13 <- 1/(1-calph^(2*diff2))
     gdd11 <- 2*(dd11^2)*diff11*(calph^(2*diff11-1))
     dd21 <- (1-calph^(2*diff2))/(1-calph^(2*diff12))
     gdd21 <- -2*(dd12^2)*(calph^(2*(diff12+diff2)-1))*(
           diff2*(calph^(-2*diff12)-1)-diff12*(calph^(-2*diff2)-1))
     gdd2 <- gdd11*dd21+dd11*gdd21
     ggdd11 <- 2*(dd11^3)*diff11*(calph^(2*diff11-2))*
            ((2*diff11-1)+(2*diff11+1)*calph^(2*diff11))
     ggdd21 <- -2*(dd12^3)*(calph^(2*(diff12+diff2)-2))*(
      diff2*(2*diff2-1)*(1+((4*diff12/(2*diff2-1))-1)*calph^(2*diff12))*(calph^(-2*diff12)-1)
      -diff12*(2*diff12-1)*(1+((4*diff12/(2*diff12-1))-1)*calph^(2*diff12))*(calph^(-2*diff2)-1))
     ggdd2 <- ggdd11*dd21+2*gdd11*gdd21+dd11*ggdd21
     IC_mat[i,j] <- dd2
     gIC_mat[i,j] <- gdd2
     ggIC_mat[i,j] <- ggdd2
    }
   }
 
    # off diagonals
     diffij <- abs(i-j)
     if (diffij==1){
      if (calph==0){
       IC_mat[i,j] <- 0
       gIC_mat[i,j] <- -1
       ggIC_mat[i,j] <- 0
      }
      else
      {
       dd3 <- -(calph^(diff))/(1-calph^(2*diff))
       dd1 <- 1/(1-calph^(2*diff))
       gdd3 <- -(dd1^2)*diff*(calph^(diff-1))*(1+calph^(2*diff))
       ggdd3 <- -(dd1^3)*diff*(calph^(diff-2))*
                 ((diff+1)*(calph^(4*diff))+6*diff*(calph^(2*diff))+(diff-1))
       IC_mat[i,j] <- dd3
       gIC_mat[i,j] <- gdd3
       ggIC_mat[i,j] <- ggdd3
     }
    }
   } # i
  } # j
  } # end ar1 correlation model


  # uniform model
  if(corstr=="uniform"){
  nn <- length(ctimes)
  C_mat <- matrix(0,nrow=nn,ncol=nn) 
  IC_mat <- C_mat; gIC_mat <- C_mat; ggIC_mat <- C_mat

  for (i in 1:length(ctimes)){
   for (j in 1:length(ctimes)){
 
   # matrix C
   if (i==j){
    C_mat[i,j] <- 1
   } else {
    C_mat[i,j] <- calph
   }

   # matrices IC, gIC and ggIC
   d_mult <- (nn-1)*(calph-1)*(calph+1/(nn-1))
   b <- calph/d_mult
   gb <- -(1+(nn-1)*(calph^2))/(d_mult^2)
   ggb <- (2*calph*(nn-1)*((nn-1)*(calph^2)+3)-2*(nn-2))/(d_mult^3)
   a <- -(1+(nn-2)*calph)/d_mult
   ga <- calph*(nn-1)*(2+(nn-2)*calph)/(d_mult^2)
   gga <- -2*(nn-1)*((calph^3)*(nn-2)*(nn-1)+3*(calph^2)*(nn-1)+1)/(d_mult^3)
   if (i==j){
    IC_mat[i,j] <- a
    gIC_mat[i,j] <- ga
    ggIC_mat[i,j] <- gga
   } else {
    IC_mat[i,j] <- b
    gIC_mat[i,j] <- gb
    ggIC_mat[i,j] <- ggb
   }
   } # i
  } # j

  } # end of uniform correlation model

  # output
  list(C_mat=C_mat,IC_mat=IC_mat,gIC_mat=gIC_mat,ggIC_mat=ggIC_mat)
 }

 # function to construct S matrix
 
 smat.work <- function(coef,categs){

  # matrix S
  categs1 <- categs-1
  S_mat <- matrix(nrow=categs1,ncol=categs1)
  for (j in 1:categs1){
   for (k in 1:categs1){
    if (j>=k){
      S_mat[j,k] <- sqrt(exp(coef[k]-coef[j]))
    }
    if (j<k){
      S_mat[j,k] <- sqrt(exp(coef[j]-coef[k]))
    }
   } # k
  } # j

  # output
  list(S_mat=S_mat)
 }

 # function to construct H matrices
   
 hmat.work <- function(mod.gee,C_mat,X_mat,S_mat,miss_data){

  # set-up structures
  Var_mat <- sqrt(mod.gee$fitted*(1-mod.gee$fitted))
  H_mat <- matrix(0,ncol=length(mod.gee$coefficients),
                    nrow=length(mod.gee$coefficients))
  gH_mat <- H_mat
  ggH_mat <- H_mat
  IS_mat <- solve(S_mat$S_mat)
  icorr_mat <- kronecker(C_mat$IC_mat,IS_mat)
  gicorr_mat <- kronecker(C_mat$gIC_mat,IS_mat)
  ggicorr_mat <- kronecker(C_mat$ggIC_mat,IS_mat)
  miss_cnt <- 0

  # loop for subjects
  for (i in 1:max(mod.gee$id)){
  
   if(sum(miss_data$miss_ind[miss_data$sub_ind==i])>0){
    # missing values
    missing <- miss_data$miss_ind[miss_data$sub_ind==i]
    miss_cnt <- miss_cnt+1
    dH_mat <- matrix(0,nrow=dim(X_mat)[2],ncol=dim(X_mat)[2])
    dgH_mat <- dH_mat
    dggH_mat <- dH_mat
   } else{

    # no missing values
    D_mat <- matrix(0,nrow=mod.gee$max.id,ncol=mod.gee$max.id)
    diag(D_mat) <- exp(mod.gee$linear.predictors[mod.gee$id==i])/
         ((1+exp(mod.gee$linear.predictors[mod.gee$id==i]))^2)
    D_mat <- t(X_mat[mod.gee$id==i,])%*%D_mat
    V_mat <- matrix(0,nrow=mod.gee$max.id,ncol=mod.gee$max.id)
    diag(V_mat) <- 1/Var_mat[mod.gee$id==i]
    Vcov_mat <- V_mat%*%icorr_mat%*%V_mat
    gVcov_mat <- V_mat%*%gicorr_mat%*%V_mat
    ggVcov_mat <- V_mat%*%ggicorr_mat%*%V_mat
    dH_mat <- D_mat%*%Vcov_mat%*%t(D_mat)
    dgH_mat <- D_mat%*%gVcov_mat%*%t(D_mat)
    dggH_mat <- D_mat%*%ggVcov_mat%*%t(D_mat)

   }

   H_mat <- H_mat+dH_mat
   gH_mat <- gH_mat+dgH_mat
   ggH_mat <- ggH_mat+dggH_mat
 
  } # for i

  # weight for missing values
  H_mat <- H_mat*max(mod.gee$id)/(max(mod.gee$id)-miss_cnt)
  gH_mat <- gH_mat*max(mod.gee$id)/(max(mod.gee$id)-miss_cnt)
  ggH_mat <- ggH_mat*max(mod.gee$id)/(max(mod.gee$id)-miss_cnt)

  # output
  list(H_mat=H_mat,gH_mat=gH_mat,ggH_mat=ggH_mat)
 }

 # function to construct G matrices

 gmat.work <- function(mod.gee,C_mat,X_mat,S_mat,categories,miss_data){

  # set-up structures
  ressers <- mod.gee$y-mod.gee$fitted
  Var_mat <- sqrt(mod.gee$fitted*(1-mod.gee$fitted))
  G_mat <- matrix(0,ncol=length(mod.gee$coefficients),
                   nrow=length(mod.gee$coefficients))
  gG_mat <- G_mat
  ggG_mat <- G_mat
  rawG_mat <- G_mat
  tind <- mod.gee$max.id/(categories-1)
  unit_mat <- matrix(0,ncol=tind,nrow=tind)
  diag(unit_mat) <- rep(1,tind)
  IS_mat <- solve(S_mat$S_mat)
  icorr_mat <- kronecker(C_mat$IC_mat,IS_mat)
  gicorr_mat <- kronecker(C_mat$gIC_mat,IS_mat)
  ggicorr_mat <- kronecker(C_mat$ggIC_mat,IS_mat)
  blockdiag <- kronecker(unit_mat,S_mat$S_mat) 
  multmat1 <- matrix(1,ncol=tind,nrow=tind)
  diag(multmat1) <- 0
  multmat2 <- matrix(1,ncol=(categories-1),nrow=(categories-1))
  multmat <- kronecker(multmat1,multmat2)
  miss_cnt <- 0

  # loop for subjects
  for (i in 1:(max(mod.gee$id))){

   if(sum(miss_data$miss_ind[miss_data$sub_ind==i])>0){
    # missing values
   
    missing <- miss_data$miss_ind[miss_data$sub_ind==i]
    miss_cnt <- miss_cnt+1
    dG_mat <- matrix(0,nrow=dim(X_mat)[2],ncol=dim(X_mat)[2])
    dgG_mat <- dG_mat
    dggG_mat <- dG_mat
    rawdG_mat <- dG_mat

   } else{
    # no missing values
   
    res_mat <- ressers[mod.gee$id==i]%*%t(ressers[mod.gee$id==i])
    D_mat <- matrix(0,nrow=mod.gee$max.id,ncol=mod.gee$max.id)
    diag(D_mat) <- exp(mod.gee$linear.predictors[mod.gee$id==i])/
         ((1+exp(mod.gee$linear.predictors[mod.gee$id==i]))^2)
    D_mat <- t(X_mat[mod.gee$id==i,])%*%D_mat
    V_mat <- matrix(0,nrow=mod.gee$max.id,ncol=mod.gee$max.id)
    diag(V_mat) <- 1/Var_mat[mod.gee$id==i]
    nblockdiag <- solve(V_mat)%*%blockdiag%*%solve(V_mat)
    raw_res_mat <- res_mat
    res_mat <- res_mat*multmat+nblockdiag
    Vcov_mat <- V_mat%*%icorr_mat%*%V_mat
    gVcov_mat <- V_mat%*%gicorr_mat%*%V_mat
    ggVcov_mat <- V_mat%*%ggicorr_mat%*%V_mat
    dG_mat <- D_mat%*%Vcov_mat%*%res_mat
    dG_mat <- dG_mat%*%Vcov_mat%*%t(D_mat)
    rawdG_mat <- D_mat%*%Vcov_mat%*%raw_res_mat
    rawdG_mat <- rawdG_mat%*%Vcov_mat%*%t(D_mat)
    dgG_mat1 <- D_mat%*%gVcov_mat%*%res_mat
    dgG_mat1 <- dgG_mat1%*%Vcov_mat%*%t(D_mat)
    dgG_mat2 <- D_mat%*%Vcov_mat%*%res_mat
    dgG_mat2 <- dgG_mat2%*%gVcov_mat%*%t(D_mat)
    dgG_mat <- dgG_mat1+dgG_mat2
    dggG_mat1 <- D_mat%*%ggVcov_mat%*%res_mat
    dggG_mat1 <- dggG_mat1%*%Vcov_mat%*%t(D_mat)
    dggG_mat2 <- D_mat%*%gVcov_mat%*%res_mat
    dggG_mat2 <- dggG_mat2%*%gVcov_mat%*%t(D_mat)
    dggG_mat3 <- D_mat%*%Vcov_mat%*%res_mat
    dggG_mat3 <- dggG_mat3%*%ggVcov_mat%*%t(D_mat)
    dggG_mat <- dggG_mat1+2*dggG_mat2+dggG_mat3
  
   }

   G_mat <- G_mat+dG_mat
   gG_mat <- gG_mat+dgG_mat
   ggG_mat <- ggG_mat+dggG_mat
   rawG_mat <- rawG_mat+rawdG_mat

  } # for i
 
  # weight for missing values
  G_mat <- G_mat*max(mod.gee$id)/(max(mod.gee$id)-miss_cnt)
  gG_mat <- gG_mat*max(mod.gee$id)/(max(mod.gee$id)-miss_cnt)
  ggG_mat <- ggG_mat*max(mod.gee$id)/(max(mod.gee$id)-miss_cnt)
  rawG_mat <- rawG_mat*max(mod.gee$id)/(max(mod.gee$id)-miss_cnt)

  # output
  list(G_mat=G_mat,gG_mat=gG_mat,ggG_mat=ggG_mat,rawG_mat=rawG_mat)
 }

 # function to update alpha

 alpha.update <- function(H_mat,G_mat){

  # eigen analysis
  heigen <- eigen(H_mat$H_mat,symmetric=TRUE)
  geigen <- eigen(G_mat$G_mat,symmetric=TRUE)

  # set-up values for sum
  gldetH_mat <- 0
  gldetG_mat <- 0
  ggldetH_mat <- 0
  ggldetG_mat <- 0

  # sum over coefficients
  for (i in 1:(dim(H_mat$H_mat)[1])){
   dheig <- t(matrix(heigen$vectors[,i],ncol=1))%*%
           H_mat$gH_mat%*%matrix(heigen$vectors[,i],ncol=1)
   gldetH_mat <- gldetH_mat+dheig/heigen$values[i]
   dgeig <- t(matrix(geigen$vectors[,i],ncol=1))%*%
           G_mat$gG_mat%*%matrix(geigen$vectors[,i],ncol=1)
   gldetG_mat <- gldetG_mat+dgeig/geigen$values[i]

   # build F matrices
   FH_mat <- matrix(0,ncol=dim(H_mat$H_mat)[1], nrow=dim(H_mat$H_mat)[1])
   FG_mat <- FH_mat
   for (k in 1:(dim(H_mat$H_mat)[1])){
     dFH_mat <- matrix(heigen$vector[,k],ncol=1)%*%
              t(matrix(heigen$vector[,k],ncol=1))
     dFG_mat <- matrix(geigen$vector[,k],ncol=1)%*%
              t(matrix(geigen$vector[,k],ncol=1))
    if (k==i){
    dFH_mat <- 0
    dFG_mat <- 0
    }
    else{
    dFH_mat <- dFH_mat/(heigen$values[k]-heigen$values[i])
    dFG_mat <- dFG_mat/(geigen$values[k]-geigen$values[i])
    }
    FH_mat <- FH_mat+dFH_mat               
    FG_mat <- FG_mat+dFG_mat 
   } # for k

  # expressions for derivatives
  ddheig <- H_mat$ggH_mat-2*H_mat$gH_mat%*%FH_mat%*%H_mat$gH_mat
  ddheig <- t(matrix(heigen$vectors[,i],ncol=1))%*%
           ddheig%*%matrix(heigen$vectors[,i],ncol=1)
  ddgeig <- G_mat$ggG_mat-2*G_mat$gG_mat%*%FG_mat%*%G_mat$gG_mat
  ddgeig <- t(matrix(geigen$vectors[,i],ncol=1))%*%
           ddgeig%*%matrix(geigen$vectors[,i],ncol=1)
  ggldetH_mat <- ggldetH_mat+(ddheig/heigen$values[i]-(dheig/heigen$values[i])^2)
  ggldetG_mat <- ggldetG_mat+(ddgeig/geigen$values[i]-(dgeig/geigen$values[i])^2) 
  } # for i 

  gVB <- gldetG_mat-2*gldetH_mat
  ggVB <- ggldetG_mat-2*ggldetH_mat
 
  # output
  list(gVB=gVB,ggVB=ggVB)
 }


 # function to update phi

 phi.update <- function(gVB,ggVB,icalph){
 
  # select transformation
  transform <- 1

  # constrain to (0,1)
  if (transform==1){
   phi <- log(icalph)-log(1-icalph)
   gphi <- exp(phi)/((1+exp(phi))^2)
   ggphi <- (exp(phi)*(1-exp(phi)))/((1+exp(phi))^3)
   gVB_phi <- gphi*gVB
   ggVB_phi <- ggphi*gVB+(gphi^2)*ggVB
   # stop phi changing too rapidly
   phi_inc <- gVB_phi/abs(ggVB_phi)
   if (abs(phi_inc)>1){
    if (phi_inc<1){phi <- phi+1}
    if (phi_inc>1){phi <- phi-1}
   }
   else{
    phi <- phi-gVB_phi/abs(ggVB_phi)
   }
   icalph <- exp(phi)/(1+exp(phi))
  }

  # constrain to (-1,1)
  if (transform==2){
   phi <- log(1+icalph)-log(1-icalph)
   gphi <- 2*exp(phi)/((1+exp(phi))^2)
   ggphi <- 2*(exp(phi)*(1-exp(phi)))/((1+exp(phi))^3)
   gVB_phi <- gphi*gVB
   ggVB_phi <- ggphi*gVB+(gphi^2)*ggVB
   # stop phi changing too rapidly
   phi_inc <- gVB_phi/abs(ggVB_phi)
   if (abs(phi_inc)>1){
    if (phi_inc<1){phi <- phi+1}
    if (phi_inc>1){phi <- phi-1}
   }
   else{
    phi <- phi-gVB_phi/abs(ggVB_phi)
   }
   icalph <- (exp(phi)-1)/(exp(phi)+1)
  }
  # output
  list(phi=phi,icalph=icalph,gVB=gVB_phi,ggVB=ggVB_phi)
 }

 # function to construct score test

 score.test <- function(mod.gee,C_mat,X_mat,S_mat,
            categories,miss_data,corstr){

  # expanded model coefficients
  excoeff <- c(mod.gee$coeff[1:(categories-1)],
      rep(mod.gee$coeff[categories:length(mod.gee$coeff)],(categories-1)))

  # set-up structures
  ressers <- mod.gee$y-mod.gee$fitted
  Var_mat <- sqrt(mod.gee$fitted*(1-mod.gee$fitted))
  W_mat <- matrix(0,ncol=length(excoeff),nrow=length(excoeff))
  U_mat <- matrix(0,ncol=1,nrow=length(excoeff))
  tind <- mod.gee$max.id/(categories-1)
  unit_mat <- matrix(0,ncol=tind,nrow=tind)
  diag(unit_mat) <- rep(1,tind)
  IS_mat <- solve(S_mat$S_mat)
  if (corstr=="independence"){
   icorr_mat <- kronecker(C_mat,IS_mat)
  }else {
   icorr_mat <- kronecker(C_mat$IC_mat,IS_mat)
  }
  blockdiag <- kronecker(unit_mat,S_mat$S_mat) 
  multmat1 <- matrix(1,ncol=tind,nrow=tind)
  diag(multmat1) <- 0
  multmat2 <- matrix(1,ncol=(categories-1),nrow=(categories-1))
  multmat <- kronecker(multmat1,multmat2)
  miss_cnt <- 0

  # loop for subjects
  for (i in 1:(max(mod.gee$id))){

   if(sum(miss_data$miss_ind[miss_data$sub_ind==i])>0){
    # missing values
    missing <- miss_data$miss_ind[miss_data$sub_ind==i]
    miss_cnt <- miss_cnt+1
    dW_mat <- matrix(0,nrow=dim(X_mat)[2],ncol=dim(X_mat)[2])
    dU_mat <- matrix(0,ncol=1,nrow=length(excoeff))
   } else{
    # no missing values
    res_mat <- ressers[mod.gee$id==i]%*%t(ressers[mod.gee$id==i])
    D_mat <- matrix(0,nrow=mod.gee$max.id,ncol=mod.gee$max.id)
    diag(D_mat) <- exp(mod.gee$linear.predictors[mod.gee$id==i])/
         ((1+exp(mod.gee$linear.predictors[mod.gee$id==i]))^2)
    D_mat <- t(X_mat[mod.gee$id==i,])%*%D_mat
    V_mat <- matrix(0,nrow=mod.gee$max.id,ncol=mod.gee$max.id)
    diag(V_mat) <- 1/Var_mat[mod.gee$id==i]
    nblockdiag <- solve(V_mat)%*%blockdiag%*%solve(V_mat)
    raw_res_mat <- res_mat
    res_mat <- res_mat*multmat+nblockdiag
    Vcov_mat <- V_mat%*%icorr_mat%*%V_mat
    dW_mat <- D_mat%*%Vcov_mat%*%res_mat
    dW_mat <- dW_mat%*%Vcov_mat%*%t(D_mat)
    dU_mat <- D_mat%*%Vcov_mat%*%ressers[mod.gee$id==i]
   }
   W_mat <- W_mat+dW_mat
   U_mat <- U_mat+dU_mat
  } # for i
 
  # output
  list(W_mat=W_mat,U_mat=U_mat)
 }

 # introductory message 
 cat(" REPOLR: 28/02/2008 version 1.1","\n")
 cat(" using gee 4.13 GEE C source version chanlib 4.12","\n")
 call <- match.call()

 # check correlation structure, maxiter, scalevalue, alpha, categories
 corstrs <- c("ar1", "uniform", "independence")
 icorstr <- as.integer(match(corstr, corstrs, -1))
 if (icorstr < 1){ 
   stop("unknown correlation structure")
 }
 cat(" working correlation model: ",corstr,"\n")
 maxiter <- as.integer(maxiter)
 scalevalue <- as.double(scalevalue)
 tol <- as.double(tol)
 alpha <- as.double(alpha)
 po.test <- as.logical(po.test)
 fixed <- as.logical(fixed)
 if (alpha<0.05 | alpha>0.95){
  stop("invalid correlation parameter (alpha)")
 }
 categories <- as.integer(categories)

 # append formula with cuts factor
 form_vars <- attr(terms.formula(as.formula(formula)),"variables")
 resp_var <- attr(terms.formula(as.formula(formula)),"response")
 term_labels <- attr(terms.formula(as.formula(formula)),"term.labels")
 if (resp_var==0) {
  stop("No response variable in formula")
 }
 resp_var <- resp_var+1
 resp_label <- as.character(form_vars[[resp_var]])
 if (is.element("cuts",term_labels)){
  stop("Term name cuts is not permitted")
 }
 formula <- as.formula(paste(resp_label,"~",
      paste(c("factor(cuts)-1",paste(term_labels,collapse="+")),collapse="+")))

 # expand data
 subjects <- as.character(subjects)
 isubject <- as.integer(match(subjects, names(as.data.frame(data)), -1))
 if (isubject < 1){ 
   stop("unknown subject name")
 }
 exdata <- ord.expand(scores=resp_label,data=data,subjects=subjects,categories=categories)

 # model matrix
 X_mat <- model.matrix(formula,data=exdata$exdata)

 # missing value indicator
 full_dat <- (categories-1)*length(times)*max(exdata$exdata$subjects)
 miss_ind <- as.vector(rep(0,full_dat),mode="numeric")
 sub_ind <- rep(1:max(exdata$exdata$subjects),each=(categories-1)*length(times))
 miss_rows <- setdiff(as.character(1:full_dat),rownames(X_mat))
 miss_ind[as.numeric(miss_rows)] <- 1
 miss_data <- as.data.frame(cbind(sub_ind,miss_ind))

 # fit gee model to get initial parameter estimates
 mod.gee <- ogee(formula,data=exdata$exdata,id=exdata$exdata$subjects,family="binomial",
     scale.fix=TRUE,scale.value=scalevalue,corstr="independence",maxiter=20,silent=TRUE)
 coeffs <- mod.gee$coeff

 # fit gee model
 iterct <- 0    # iteration count
 ogstop <- 1    # stopping criterion

 while (ogstop>tol & iterct<=maxiter){

  # update
  if (corstr=="independence"){
   C_mat <- matrix(0,nrow=length(times),ncol=length(times)) 
   diag(C_mat) <- rep(1,length(times))
  } else{
   C_mat <- cmat.work(calph=alpha,ctimes=times,corstr=corstr) 
  }
  categories1 <- categories-1
  S_mat <- smat.work(coeffs[1:categories1],categories)
  if (corstr=="independence"){
   R_mat <- kronecker(C_mat,S_mat$S_mat)
  } else {
   R_mat <- kronecker(C_mat$C_mat,S_mat$S_mat)
  }

  if (corstr!="independence"){
  
  if (fixed==FALSE){
   # h matrices
   H_mat <- hmat.work(mod.gee=mod.gee,C_mat=C_mat,X_mat=X_mat,S_mat=S_mat,miss_data=miss_data)
   # g matrices
   G_mat <- gmat.work(mod.gee=mod.gee,C_mat=C_mat,X_mat=X_mat,S_mat=S_mat,
          categories=categories,miss_data=miss_data)
   # eigen vectors and values
   aupdate <- alpha.update(H_mat=H_mat,G_mat=G_mat)
   pupdate <- phi.update(gVB=aupdate$gVB,ggVB=aupdate$ggVB,icalph=alpha)

   # update alpha
   alpha <- pupdate$icalph
   if (alpha>0.95){alpha <- 0.95}
   if (alpha<0.05){alpha <- 0.05}
  }
  if (fixed==TRUE){
   alpha <- alpha
  }
  }

  # fit GEE
  mod.gee <- ogee(formula,data=exdata$exdata,id=exdata$exdata$subjects,
     family="binomial",b=as.numeric(coeffs),scale.fix=TRUE,scale.value=scalevalue,
     corstr="fixed",maxiter=10,R=R_mat,silent=TRUE)

  # stopping criterion
  crit <- abs(1-sqrt(sum((coeffs/mod.gee$coeff)^2)/length(coeffs)))
  coeffs <- mod.gee$coeff

  # monitor
  iterct <- iterct+1
  if (iterct==1){
    cat(format("iter",digits=1,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=8),
        "\t",format("crit",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
        "\t",format("alpha",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
        "\t",format("grad1",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
        "\t",format("grad2",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),"\n")
  }

  if (corstr=="independence"){
     cat(format(round(iterct,0),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=8),
     "\t",format(round(crit,5),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format("-",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format("-",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format("-",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),"\n")
  } else {
    if (fixed==FALSE){
    cat(format(round(iterct,0),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=8),
     "\t",format(round(crit,5),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format(round(alpha,5),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format(round(pupdate$gVB,5),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format(round(pupdate$ggVB,5),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),"\n")
    }
    if (fixed==TRUE){
     cat(format(round(iterct,0),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=0,width=8),
     "\t",format(round(crit,5),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format(round(alpha,5),digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format("-",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),
     "\t",format("-",digits=5,trim=TRUE,justify="right",scientific=FALSE,nsmall=5,width=8),"\n")
    }
  }
  if (iterct<=4){
   ogstop <- 1
  }else{
   ogstop <- crit
  }

  } # end of while


  # score test
  if (po.test==TRUE){
   exformula <- as.formula(paste(resp_label,"~",
            "factor(cuts)*(",paste(term_labels,collapse="+"),
                  paste(")-1",sep=""),sep=""))
   exX_mat <- model.matrix(exformula,data=exdata$exdata)
   s.test <- score.test(mod.gee=mod.gee,C_mat=C_mat,X_mat=exX_mat,S_mat=S_mat,
          categories=categories,miss_data=miss_data,corstr=corstr)
   test_stat <- t(s.test$U_mat)%*%solve(s.test$W_mat)%*%s.test$U_mat
   test_df <- (categories-2)*(length(mod.gee$coeff)-(categories-1))
   test_chi <- pchisq(test_stat,test_df,lower.tail=FALSE)
   cat("\n")
   cat(" Score test = ",round(test_stat,3),": d.f. ",round(test_df,0),
        " and p-value ",round(test_chi,3),"\n")
  }

  # present CI on original scale
  if (corstr!="independence" & fixed==FALSE){
  if (pupdate$ggVB>0){
   phi <- log(alpha)-log(1-alpha)
   phi_u <- phi+1.96*(1/sqrt(pupdate$ggVB))
   phi_l <- phi-1.96*(1/sqrt(pupdate$ggVB))
   alpha_u <- exp(phi_u)/(1+exp(phi_u))
   alpha_l <- exp(phi_l)/(1+exp(phi_l))
   if (po.test!=TRUE){cat("\n")}
   cat(" Correlation parameter (95% CI):",round(alpha,3)," (",round(alpha_l,3),
                     ", ",round(alpha_u,3),")","\n")
  } else {
   warning(" grad2 < 0: minimum for alpha not achieved")
  }
 }

 # fitted model
 mod.corr <- list()
 mod.corr$subjects <- subjects
 mod.corr$categories <- categories
 mod.corr$times <- times
 mod.corr$corstr <- corstr
 mod.corr$crit <- crit
 mod.corr$score.test <- test_stat
 mod.corr$iter <- iterct
 mod.corr$scale <- scalevalue
 mod.corr$data <- call[["data"]]
 if (corstr!="independence"){
  mod.corr$alpha <- alpha
  mod.corr$grad1 <- pupdate$gVB
  mod.corr$grad2 <- pupdate$ggVB
 }
 list(gee=mod.gee,corr=mod.corr)
}

