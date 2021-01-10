library(fGarch)

start.date <- '2009-01-01'
prc.history <- read.csv('data/^DJI.csv', header = TRUE)
prices <- prc.history[prc.history['Date'] >= start.date, 'Adj.Close']
prc.rtn <- log(prices[-1]) - log(prices[-length(prices)])
gm <- garchFit(~garch(1,1), data = prc.rtn, include.mean = FALSE)

res <- residuals(gm, standardize = TRUE)
res[130] * gm@sigma.t[130]
# model residuals are not normalized.
gm@residuals[130]
prc.rtn[130]

# Verify GARCH
gm@fit$par[1] + gm@fit$par[2] * gm@residuals[129]^2 + gm@fit$par[3] * gm@sigma.t[129]^2
gm@sigma.t[130]^2
