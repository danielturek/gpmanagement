

library(nimble)

gpPred <- nimble::nimbleFunction(
    setup = function(model, params) {
        calcNodes <- model$getDependencies(params, determOnly = TRUE)
        nPred <- dim(model$SigPP)[1]
        E <- array(0, c(nPred, 1))
        C <- array(0, c(nPred, nPred))
    },
    run = function(samples = double(2)) {
        E <<- E * 0
        C <<- C * 0
        nSamples <- dim(samples)[1]
        for(i in 1:nSamples) {
            values(model, params) <<- samples[i,]
            calculate(model, calcNodes)
            intermediate <- model$SigPO %*% inverse(model$SigOO)
            Etemp <- intermediate %*% asCol(model$yObs)
            Ctemp <- model$SigPP - intermediate %*% t(model$SigPO)
            E <<- E + Etemp
            C <<- C + Ctemp
        }
        E <<- E / nSamples
        C <<- C / nSamples
    },
    methods = list(
        getE = function() { returnType(double(1)); return(E[,1]) },
        getC = function() { returnType(double(2)); return(C[, ]) }
    )
)

