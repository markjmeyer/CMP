#### NICOMPR Simulation ####
#   by:         Meyer, MJ, Graye, A, & Sellers, KF
#   modified:   04/06/23
#
#   Underdispersed CMP, lambda = 3, nu = 2

#### load code ####
source('nicompr.R')
options(mc.cores = parallel::detectCores())
library(COMPoissonReg)

#### settings ####
B             <- 100
n             <- c(25, 75, 125)
lambda        <- 3
nu            <- 2
x             <- c(2, 0)
priors_c1     <- list(a = 1, b = 1, c = 1)        # add one
priors_c2     <- list(a = sum(x), b = sum(log(factorial(x))), c = length(x)) # add two subjects, one with x_i = large, x_i = 0
priors_c01    <- list(a = 0.1, b = 0.1, c = 0.1)
priors_c001   <- list(a = 0.01, b = 0.01, c = 0.01)

#### storage arrays ####
loos      <- array(NA, dim = c(B, 6, length(n)))
waics     <- array(NA, dim = c(B, 6, length(n)))
lamMat    <- array(NA, dim = c(B, 6, length(n)))
nuMat     <- array(NA, dim = c(B, 6, length(n)))
lamVar    <- array(NA, dim = c(B, 6, length(n)))
nuVar     <- array(NA, dim = c(B, 6, length(n)))
cov_lam   <- array(NA, dim = c(B, 6, length(n)))
cov_nu    <- array(NA, dim = c(B, 6, length(n)))
diverge   <- array(NA, dim = c(B, 6, length(n)))
rhat_lam  <- array(NA, dim = c(B, 6, length(n)))
rhat_nu   <- array(NA, dim = c(B, 6, length(n)))

#### simulation run ####
for(i in 1:length(n)){
  for(b in 1:B){

    # i <- b <- 2
    # sts <- proc.time()
    
    ##### generate data #####
    set.seed(b)
    X   <- rcmp(n[i], lambda = lambda, nu = nu)
    
    ##### conjugate models #####
    mod_c1    <- suppressWarnings(cmp(X, priors = priors_c1,
                                      type = 'conjugate'))
    mod_c2    <- suppressWarnings(cmp(X, priors = priors_c2,
                                      type = 'conjugate'))
    mod_c01   <- suppressWarnings(cmp(X, priors = priors_c01,
                                      type = 'conjugate'))
    mod_c001  <- suppressWarnings(cmp(X, priors = priors_c001,
                                      type = 'conjugate'))
    
    ##### flat model #####
    mod_f     <- suppressWarnings(cmp(X, priors = list(a = 0, b = 0, c = 0),
                                      type = 'flat'))
    
    ##### Jeffreys' model #####
    mod_j     <- suppressWarnings(cmp(X, type = 'Jeffreys',
                                      control = list(adapt_delta = 0.999,
                                                     stepsize = 0.01,
                                                     max_treedepth = 15)))
    
    loos[b,1,i]       <- loo(mod_c1)
    loos[b,2,i]       <- loo(mod_c2)
    loos[b,3,i]       <- loo(mod_c01)
    loos[b,4,i]       <- loo(mod_c001)
    loos[b,5,i]       <- loo(mod_f)
    loos[b,6,i]       <- loo(mod_j)
    
    waics[b,1,i]      <- WAIC(mod_c1)
    waics[b,2,i]      <- WAIC(mod_c2)
    waics[b,3,i]      <- WAIC(mod_c01)
    waics[b,4,i]      <- WAIC(mod_c001)
    waics[b,5,i]      <- WAIC(mod_f)
    waics[b,6,i]      <- WAIC(mod_j)
    
    lamMat[b,1,i]     <- coef(mod_c1)[1] 
    lamMat[b,2,i]     <- coef(mod_c2)[1] 
    lamMat[b,3,i]     <- coef(mod_c01)[1] 
    lamMat[b,4,i]     <- coef(mod_c001)[1] 
    lamMat[b,5,i]     <- coef(mod_f)[1] 
    lamMat[b,6,i]     <- coef(mod_j)[1] 
    
    nuMat[b,1,i]      <- coef(mod_c1)[2] 
    nuMat[b,2,i]      <- coef(mod_c2)[2] 
    nuMat[b,3,i]      <- coef(mod_c01)[2] 
    nuMat[b,4,i]      <- coef(mod_c001)[2] 
    nuMat[b,5,i]      <- coef(mod_f)[2] 
    nuMat[b,6,i]      <- coef(mod_j)[2] 
    
    Post_c1        <- posts(mod_c1, pars = c('lambda', 'nu'))
    Post_c2        <- posts(mod_c2, pars = c('lambda', 'nu'))
    Post_c01       <- posts(mod_c01, pars = c('lambda', 'nu'))
    Post_c001      <- posts(mod_c001, pars = c('lambda', 'nu'))
    Post_f         <- posts(mod_f, pars = c('lambda', 'nu'))
    Post_j         <- posts(mod_j, pars = c('lambda', 'nu'))
    
    lamVar[b,1,i]     <- var(Post_c1$lambda)
    lamVar[b,2,i]     <- var(Post_c2$lambda)
    lamVar[b,3,i]     <- var(Post_c01$lambda)
    lamVar[b,4,i]     <- var(Post_c001$lambda)
    lamVar[b,5,i]     <- var(Post_f$lambda)
    lamVar[b,6,i]     <- var(Post_j$lambda)
    
    nuVar[b,1,i]      <- var(Post_c1$nu)
    nuVar[b,2,i]      <- var(Post_c2$nu)
    nuVar[b,3,i]      <- var(Post_c01$nu)
    nuVar[b,4,i]      <- var(Post_c001$nu)
    nuVar[b,5,i]      <- var(Post_f$nu)
    nuVar[b,6,i]      <- var(Post_j$nu)
    
    ci_c1           <- credint(mod_c1)
    ci_c2           <- credint(mod_c2)
    ci_c01          <- credint(mod_c01)
    ci_c001         <- credint(mod_c001)
    ci_f            <- credint(mod_f)
    ci_j            <- credint(mod_j)
    
    cov_lam[b,1,i]    <- 1*(ci_c1$lambda[1] < lambda & ci_c1$lambda[2] > lambda)
    cov_lam[b,2,i]    <- 1*(ci_c2$lambda[1] < lambda & ci_c2$lambda[2] > lambda)
    cov_lam[b,3,i]    <- 1*(ci_c01$lambda[1] < lambda & ci_c01$lambda[2] > lambda)
    cov_lam[b,4,i]    <- 1*(ci_c001$lambda[1] < lambda & ci_c001$lambda[2] > lambda)
    cov_lam[b,5,i]    <- 1*(ci_f$lambda[1] < lambda & ci_f$lambda[2] > lambda)
    cov_lam[b,6,i]    <- 1*(ci_j$lambda[1] < lambda & ci_j$lambda[2] > lambda)
    
    cov_nu[b,1,i]     <- 1*(ci_c1$nu[1] < nu & ci_c1$nu[2] > nu)
    cov_nu[b,2,i]     <- 1*(ci_c2$nu[1] < nu & ci_c2$nu[2] > nu)
    cov_nu[b,3,i]     <- 1*(ci_c01$nu[1] < nu & ci_c01$nu[2] > nu)
    cov_nu[b,4,i]     <- 1*(ci_c001$nu[1] < nu & ci_c001$nu[2] > nu)
    cov_nu[b,5,i]     <- 1*(ci_f$nu[1] < nu & ci_f$nu[2] > nu)
    cov_nu[b,6,i]     <- 1*(ci_j$nu[1] < nu & ci_j$nu[2] > nu)
    
    diverge[b,1,i]    <- mod_c1$hmc$numDiverge
    diverge[b,2,i]    <- mod_c2$hmc$numDiverge
    diverge[b,3,i]    <- mod_c01$hmc$numDiverge
    diverge[b,4,i]    <- mod_c001$hmc$numDiverge
    diverge[b,5,i]    <- mod_f$hmc$numDiverge
    diverge[b,6,i]    <- mod_j$hmc$numDiverge
    
    rhat_lam[b,1,i]   <- mod_c1$diagnostics$Rhat$lambda
    rhat_lam[b,2,i]   <- mod_c2$diagnostics$Rhat$lambda
    rhat_lam[b,3,i]   <- mod_c01$diagnostics$Rhat$lambda
    rhat_lam[b,4,i]   <- mod_c001$diagnostics$Rhat$lambda
    rhat_lam[b,5,i]   <- mod_f$diagnostics$Rhat$lambda
    rhat_lam[b,6,i]   <- mod_j$diagnostics$Rhat$lambda
    
    rhat_nu[b,1,i]   <- mod_c1$diagnostics$Rhat$nu
    rhat_nu[b,2,i]   <- mod_c2$diagnostics$Rhat$nu
    rhat_nu[b,3,i]   <- mod_c01$diagnostics$Rhat$nu
    rhat_nu[b,4,i]   <- mod_c001$diagnostics$Rhat$nu
    rhat_nu[b,5,i]   <- mod_f$diagnostics$Rhat$nu
    rhat_nu[b,6,i]   <- mod_j$diagnostics$Rhat$nu
    
    # ets <- proc.time() - sts
    
  }
}


#### save output ####
save.image('udis_sim.RData')

#### tabular evals ####

modNames  <- c('conj1', 'conj2', 'conj01', 'conj001', 'flat', 'jeff')

##### bias #####
matrix(apply(lamMat - lambda, c(2,3), median), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))
matrix(apply(nuMat - nu, c(2,3), median), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))

##### MSE #####
matrix(apply((lamMat - lambda)^2, c(2,3), median), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))
matrix(apply((nuMat - nu)^2, c(2,3), median), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))

##### variance #####
matrix(apply(lamVar, c(2,3), median), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))
matrix(apply(nuVar, c(2,3), median), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))

##### coverage #####
matrix(apply(cov_lam, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))
matrix(apply(cov_nu, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))

##### divergent chains ####
matrix(apply(diverge, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))
matrix(apply(diverge, c(2,3), mean)/4000, nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))


##### Rhat #####
matrix(apply(rhat_lam, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))
matrix(apply(rhat_nu, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))

matrix(apply(rhat_lam < 1.01, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))
matrix(apply(rhat_nu < 1.01, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))


##### WAICs #####
matrix(apply(waics, c(2,3), mean), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))


##### loos #####
matrix(apply(loos, c(2,3), median), nrow = 6, byrow = FALSE,
       dimnames = list(Models = modNames, `Sample Sizes` = n))

