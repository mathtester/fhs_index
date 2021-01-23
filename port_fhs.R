library(fGarch)
library(MASS)
options(warn=-1)

# The total run time is about 20 min on a PC.

calc.var <- function(price.df, garch.period, corr.period,
                     var.period, var.rank, idx.cnt) {
  # price.df is a data frame with dates as row names.
  # Each column is a series of historical prices.
  # idx.cnt is the number of indices in the columns of price.df.
  # It is assumed the stock indices are in the first columns.
  symbols <- colnames(price.df)
  stocks <- symbols[(idx.cnt + 1):length(symbols)]  # individual stock tickers.
  rtn.df <- price.df[-1, ] / price.df[-nrow(price.df), ] - 1
  garch.cols <- c('omega', 'alpha', 'beta', 'sigma', 'residual', 'return')
  subcols <- c(garch.cols, 'var')
  sym.cols <- outer(symbols, subcols, FUN = 'paste', sep = '.')
  sym.cols <- t(sym.cols)
  dim(sym.cols) <- NULL    # Flatten to a 1D vector.
  all.cols <- c('date', 'corr.th', 'corr.ex', 'port.var', sym.cols)
  garch.df <- data.frame(matrix(nrow = 0, ncol = length(all.cols)))
  colnames(garch.df) <- all.cols
  garch.days <- nrow(rtn.df) - garch.period + 1
  for (j in 1:garch.days) {
    selected = rtn.df[j:(j+garch.period-1), ]
    curr.date <- rownames(rtn.df)[j+garch.period-1]
    garch.df[j, 'date'] <- curr.date
    for (s in symbols) {
      # skip outputs of garchFit().
      capture.output( {
        gm <- garchFit(~garch(1,1), data = selected[, s], include.mean = FALSE)
      } )
      # This is the forecast sigma using data up to prior period.
      sig <- gm@sigma.t[garch.period]
      # Residuals are not normalized.
      rtn <- gm@residuals[garch.period]
      garch.df[j, paste(s, garch.cols, sep = '.')] <- 
        c(gm@fit$par, sig, rtn / sig, rtn)
    }
    
    if (j > var.period) {
      # Now compute average correlations.
      norm.rtn <- as.matrix(garch.df[(j-corr.period+1):j,
                                     paste(stocks, 'residual', sep = '.')],
                            nrow = corr.period)
      corr.mat <- cor(norm.rtn)
      d <- nrow(corr.mat)
      garch.df[j, 'corr.ex'] <- (sum(corr.mat) - d) / d / (d-1)
      last.rtn <- matrix(norm.rtn[nrow(norm.rtn), ], nrow = 1)
      garch.df[j, 'corr.th'] <- (sum(t(last.rtn) %*% last.rtn) - sum(last.rtn^2))/d/(d-1)
      
      # Compute FHS VaR.
      residual <- garch.df[(j-var.period):(j-1), paste(symbols, 'residual', sep = '.')]
      sigma <- matrix(as.numeric(rep(garch.df[j, paste(symbols, 'sigma', sep = '.')],
                                 var.period)), nrow = var.period, byrow = TRUE)
      pl.df <- residual * sigma
      sorted.pls <- apply(pl.df, 2, 'sort', partial = var.rank)
      garch.df[j, paste(symbols, 'var', sep = '.')] <- -sorted.pls[var.rank, ]
      
      # Compute portfolio VaR.
      prices <- data.matrix(price.df[curr.date, stocks])
      port.pl <- data.matrix(pl.df[, (idx.cnt+1):length(symbols)]) %*%
        t(prices) / sum(prices)
      port.pl <- sort(port.pl, partial = var.rank)
      garch.df[j, 'port.var'] <- -port.pl[var.rank]
    }
  }
  return(garch.df)
}


load.dji <- function() {
  # This file is obtained from DJI component at the end of 2020,
  # with DOW removed due to its short history.
  compo <- read.csv('input/components.csv', header = TRUE)
  start.date <- '2018-01-01'
  tickers <- compo[, 1]
  price.df <- NULL
  for (t in tickers) {
    prc.history <- read.csv(paste('input/', t, '.csv', sep=''), header = TRUE)
    if (is.null(price.df)) {
      price.df <- as.data.frame(
        prc.history[prc.history['Date'] >= start.date, 'Adj.Close'],
        row.names = prc.history[prc.history['Date'] >= start.date, 'Date'])
    } else {
      price.df <- cbind(price.df, prc.history[prc.history['Date'] >= start.date, 'Adj.Close'])
    }
  }
  colnames(price.df) <- tickers
  
  port.cor <- cor(price.df)
  d <- nrow(port.cor)
  avg.cor <- (sum(port.cor) - d) / d / (d-1)
  print(sprintf('%s Average return correlation in DJI is %f', Sys.time(), avg.cor))
  
  price.df['index'] <- rowSums(price.df)
  prc.history <- read.csv('input/^DJI.csv', header = TRUE)
  price.df['DJI'] <- prc.history[prc.history['Date'] >= start.date, 'Adj.Close']
  price.df <- price.df[c(d+2, d+1, 1:d)]
  
  return(price.df)
}


simu.port <- function() {
  nstock <- 30
  regimes <- c( 800,   50,   50,   50,   50,   50,   50,  50,  100 )
  sigmas <-  c(0.01, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02 )
  correls <- c( 0.3,  0.3,  0.9,  0.9,  0.3,    0,    0,  0.9, 0.4 )
  init.prc <- rep(100, nstock)
  price.df <- data.frame(matrix(init.prc, nrow = 1))
  parm.df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(parm.df) <- c('sigma', 'correl')
  ones <- matrix(rep(1, nstock^2), nrow = nstock)
  mu <- rep(0, nstock)
  cnt <- 1
  for (j in 1:length(regimes)) {
    corr.mat <- correls[j] * ones + (1 - correls[j]) * diag(nstock)
    rtns <- mvrnorm(regimes[j], mu, corr.mat)
    for (d in 1:regimes[j]) {
      parm.df[cnt, ] <- c(sigmas[j], correls[j])
      cnt <- cnt + 1
      price.df[cnt, ] <- price.df[cnt-1, ] * (1 + rtns[d, ] * sigmas[j])
    }
  }
  price.df['index'] <- rowSums(price.df)
  price.df <- price.df[c(nstock+1, 1:nstock)]
  return(list(price.df, parm.df))
}


print(paste(Sys.time(), 'Start'))

# Construct price data frame for all symbols and the portfolio.
garch.period <- 250
corr.period <- 10
set.seed(5)

if (TRUE) {
  # Use historical DJI component prices.
  var.period <- 250
  var.rank <- 2
  price.df <- load.dji()
  fhs <- calc.var(price.df, garch.period, corr.period, var.period, var.rank, 2)
  outfile <- 'output/fhs_dji.csv'
  write.csv(fhs, file = outfile, quote = FALSE, sep = ',')
  print(paste(Sys.time(), 'Finished DJI. See', outfile))
}


if (TRUE) {
  # Use simulated prices.
  var.period <- 500
  var.rank <- 5
  x <- simu.port()
  price.df <- x[[1]]
  parm.df <- x[[2]]
  fhs <- calc.var(price.df, garch.period, corr.period, var.period, var.rank, 1)
  fhs <- cbind(parm.df[garch.period:nrow(parm.df), ], fhs)
  fhs <- fhs[, c(3, 1, 2, 4:ncol(price.df))]
  outfile <- 'output/fhs_sim.csv'
  write.csv(fhs, file = outfile, quote = FALSE, row.names = FALSE, sep = ',')
  print(paste(Sys.time(), 'Finished simulation. See', outfile))
}

