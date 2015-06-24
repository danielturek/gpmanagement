testthat::test_that("A basic GP example from Daniel", {

  library("nimble")
  set.seed(0)
  x <- 1:100
  y <- sin(x/5) + rnorm(100, 0.1)
  ind <- sort(sample(1:100, 40))
  xObs <- x[ind]                ## input
  yObs <- y[ind]                ## input
  xPred <- c(x, 101:120)        ## input

  fit <- gp_setup(xObs, yObs, xPred)

  Cmcmc <- fit$Cmcmc 
  Cpred <- fit$Cpred
  Cmodel <- fit$Cmodel
  
  system.time(Cmcmc$run(100000))

  ##About 25 seconds on my computer
  ##That can be improved if necessary

  samples <- as.matrix(Cmcmc$mvSamples)
  
  ## basic sanity check
  testthat::expect_identical(Cmodel$getNodeNames(topOnly = TRUE), colnames(samples))

  #### predict from GP model using posterior MCMC samples
  system.time(Cpred$run(samples))

  ## About 40 seconds on my computer
  ## Again, could streamline the gpPred() function if necessary

  #### extract predictions: E and C
  E <- Cpred$getE()
  C <- Cpred$getC()
  
  if(interactive()){
    plot(xObs, yObs, type='b', pch=19, xlim=range(c(xObs,xPred)), ylim=range(c(yObs,E)))
    points(xPred, E, pch=20, col='red')
    segments(x0=xPred, y0=E-sqrt(diag(C)), y1=E+sqrt(diag(C)), col='red')
  }
  
  testthat::expect_is(E, "numeric")
  testthat::expect_is(C, "matrix")
})
