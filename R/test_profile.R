fs = list.files()
for (f in fs) if (regexpr("test",f)[1]!=1) source(f)

library(Rcpp)
library(RcppArmadillo)
library(DiceKriging)
library(KrigInv)
library(grid)
library(mnormt)

sourceCpp("../src/fast_suroptinv.cpp")
sourceCpp("../src/RcppExports.cpp")


library(KrigInv)
myfun <- branin_robinv
d <- 3

set.seed(8)

n0 <- 30
T <- 10
opt.index <- c(3)
inv.index <- c(1,2)
lower <- rep(0,times=d)
upper <- rep(1,times=d)

design <- matrix(runif(d*n0),nrow=n0)
response <- myfun(design)
model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2")

integcontrol <- list(distrib = "surnew",n.points = 20,finaldistrib="surnew",
                     n.candidates=50,nsimu=1000,n.optpoints = 50,
                     choose_optpoints=TRUE,n.optpoints.candidates=500)

obj <- integration_design_robinv(integcontrol = integcontrol,d=d,lower=lower,upper=upper,
                                 opt.index=opt.index,inv.index=inv.index,model=model,T=T)

# one try with discrete optimization:
optimcontrol <- list(method="discrete")

result <- max_surnew_robinv(lower = lower,upper = upper,optimcontrol = optimcontrol,
                            opt.index = opt.index,inv.index = inv.index,
                            integration.param = obj,T = T,model = model)

result$par
result$value

# one try with genoud optimization:
optimcontrol <- list(method="genoud",pop.size = 200,max.generations=3)

result2 <- max_surnew_robinv(lower = lower,upper = upper,optimcontrol = optimcontrol,
                             opt.index = opt.index,inv.index = inv.index,
                             integration.param = obj,T = T,model = model)

result2$par
result2$value
