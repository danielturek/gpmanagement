
#' gp
#' 
#' gp 
#' @param xObs observed x values
#' @param yObs observed y values
#' @param xPred predicted x values
#' @return a list object containing handles to the compiled MCMC routine, prediction routine, and nimble model object
#' @import nimble
#' @export
gp_setup <- function(xObs, yObs, xPred){

##Some initial processing
  nObs <- length(xObs)
  nPred <- length(xPred)
  f <- function(xi, xj) (xi-xj)^2
  diffOO <- outer(xObs,  xObs,  f)
  diffPP <- outer(xPred, xPred, f)
  diffPO <- outer(xPred, xObs,  f)

##### GP model
##GP model defined here.  Basically the same as your original, but I renamed some paramters more to my liking =)
##This is *not* as elegant as I had hoped.  It still requires repeating (essentially) the same code three times.  I discovered some limitations of NIMBLE while trying other approaches.

##Bottom line:  
##- This achieves relatively nice simplicity.  
##- There's an efficiency hit to the MCMC sampling, but it doesn't seem too bad.  
##- This works.  

  code <- nimbleCode({
     rho ~ dgamma(10, 1)
     sigGP ~ dunif(0, 1e5)
     sigOE ~ dunif(0, 1e5)
     SigOO[1:nObs, 1:nObs ] <- sigGP^2*exp(-1/2*diffOO[1:nObs, 1:nObs ]/rho^2) + sigOE^2*IOO[1:nObs, 1:nObs ]
     SigPP[1:nPred,1:nPred] <- sigGP^2*exp(-1/2*diffPP[1:nPred,1:nPred]/rho^2) + sigOE^2*IPP[1:nPred,1:nPred]
     SigPO[1:nPred,1:nObs ] <- sigGP^2*exp(-1/2*diffPO[1:nPred,1:nObs ]/rho^2)
     yObs[1:nObs] ~ dmnorm(mu[1:nObs], cov = SigOO[1:nObs,1:nObs])
  })

  constants <- list(nObs=nObs, nPred=nPred, diffOO=diffOO, diffPP=diffPP, diffPO=diffPO,
                    IOO=diag(nObs), IPP=diag(nPred), mu=rep(0,nObs))

  data <- list(yObs = yObs)
  inits <- list(rho = 1, sigGP = 1, sigOE = 1)
  Rmodel <- nimbleModel(code, constants, data, inits)

  ##### Create, compile, and run

  #### MCMC specification with no samplers
  spec <- configureMCMC(Rmodel, nodes = NULL)

  #### this will be used for some checking, and setting up block sampler:
  params <- Rmodel$getNodeNames(topOnly = TRUE)

  #### NOTE: the next line shouldn't work for you, since the MCMC API has changed.
  #### you can rebuild from source off nimble branch 'devel'.
  #### otherwise, using last public release of NIMBLE (0.3-1), use this line instead:
  #### spec$addSampler('RW_block', control = list(targetNodes = params))
  spec$addSampler(params, 'RW_block')

  ##We can debate about univariate vs. block sampling at some point.

  #### MCMC function
  Rmcmc <- buildMCMC(spec)

  #### GP prediction function
  #### also uses the 'params' variable for specialization
  Rpred <- gpPred(Rmodel, params)

  #### compile everything
  Cmodel <- compileNimble(Rmodel)
  Cmcmc  <- compileNimble(Rmcmc, project = Rmodel)
  Cpred  <- compileNimble(Rpred, project = Rmodel)

  fit <- list(Cmcmc=Cmcmc, Cpred=Cpred, Cmodel=Cmodel)
  fit
}


#' gp_mcmc
#'
#' gp_mcmc
#'
#' @param Cmcmc A compiled NIMBLE MCMC function object
#' @export
gp_mcmc <- function(Cmcmc){

  ## Consider as a separate call?
  #### MCMC sampling
  #### 100,000 iterations
  system.time(Cmcmc$run(100000))

  Cmcmc
}


