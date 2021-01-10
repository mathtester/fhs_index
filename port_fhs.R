library(fGarch)
options(warn=-1)

print(paste(Sys.time(), 'Start'))


calc.var <- function(price.df, garch.period, corr.period,
                     var.period, var.rank, idx.cnt) {
  # price.df is a data frame with dates as row names.
  # Each column is a series of historical prices.
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
      # garchFit dump un
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


# Construct price data frame for all symbols and the portfolio.
start.date <- '2018-01-01'
garch.period <- 250
corr.period <- 5
var.period <- 250
var.rank <- 2
idx.cnt <- 2
# This file is obtained from DJI component at the end of 2020,
# with DOW removed due to its short history.
compo <- read.csv('input/components.csv', header = TRUE)
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
print(sprintf('%s Average return correlation is %f', Sys.time(), avg.cor))

price.df['index'] <- rowSums(price.df)
prc.history <- read.csv('input/^DJI.csv', header = TRUE)
price.df['DJI'] <- prc.history[prc.history['Date'] >= start.date, 'Adj.Close']
price.df <- price.df[c(d+2, d+1, 1:d)]
fhs <- calc.var(price.df, garch.period, corr.period, var.period, var.rank, idx.cnt)
outfile <- 'output/fhs.csv'
write.csv(fhs, file = outfile, quote = FALSE, sep = ',')


print(paste(Sys.time(), 'End. See', outfile))
