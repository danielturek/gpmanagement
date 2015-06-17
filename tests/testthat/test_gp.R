test_that("A basic GP example from Daniel", {

  x <- 1:100
  y <- sin(x/5) + rnorm(100, 0.1)
  ind <- sort(sample(1:100, 40))
  xObs <- x[ind]                ## input
  yObs <- y[ind]                ## input
  xPred <- c(x, 101:120)        ## input

  fit <- gp(xObs, yObs, xPred)

  Cmcmc <- fit$Cmcmc 
  Cpred <- fit$Cpred

  system.time(Cmcmc$run(100000))

  ##About 25 seconds on my computer
  ##That can be improved if necessary

  samples <- as.matrix(Cmcmc$mvSamples)
  expect_identical(params, colnames(samples))

  #### predict from GP model using posterior MCMC samples
  system.time(Cpred$run(samples))

  ## About 40 seconds on my computer
  ## Again, could streamline the gpPred() function if necessary

  #### extract predictions: E and C
  E <- Cpred$getE()
  C <- Cpred$getC()

  expect_is(NULL, "NULL")
})
