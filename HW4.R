# Homework 4

# load the data
load('data/comp.RData')
load('data/smr.Rdata')

# calculate pyear for accounting data
comp$pyear <- comp$year + 1
comp$year <- NULL
comp$month <- NULL
comp$mktcap <- NULL

# prepare the returns data
crsp.clean$pmonth <- ifelse(crsp.clean$month < 7, crsp.clean$month + 6, crsp.clean$month - 6)
returns <- crsp.clean[, c('PERMNO', 'pyear', 'pmonth', 'mktcap', 'lagmktcap', 'momentum', 'reversal', 'RET', 'EXCHCD')]
names(returns) <- tolower(names(returns))

# merge the accounting and returns data
compcrsp <- merge(comp, returns, by=c('permno', 'pyear'), all.y=TRUE)
compcrsp <- compcrsp[compcrsp$pyear %in% comp$pyear,]
compcrsp <- compcrsp[order(compcrsp$permno, compcrsp$pyear, compcrsp$pmonth),]

# accounting variables
vars <- c('btm', 'roa', 'agr', 'nsi', 'acc', 'size', 'mom', 'rev')
annual <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)

# rename size, mom, rev
compcrsp$size <- compcrsp$mktcap
compcrsp$mom <- compcrsp$momentum
compcrsp$rev <- compcrsp$reversal
compcrsp$mktcap <- NULL
compcrsp$momentum <- NULL
compcrsp$reversal <- NULL

# take log of net stock issues
compcrsp$nsi <- log(compcrsp$nsi)

calc.breakpoints <- function(breakpoints, compcrsp, var, probabilities, annual) {
    if (annual) {
        rows <- compcrsp$exchcd==1 & compcrsp$pyear==breakpoints[1] & compcrsp$pmonth==1
        return(quantile(compcrsp[rows, var], probabilities, na.rm=TRUE))
    } else {
        rows <- compcrsp$exchcd==1 & compcrsp$pyear==breakpoints[1] & compcrsp$pmonth==breakpoints[2]
        return(quantile(compcrsp[rows, var], probabilities, na.rm=TRUE))
    }
}

# calculate breakpoints for BTM, ROA, asset growth, net stock issues, accruals, size, momentum and reversals
compcrsp <- compcrsp[compcrsp$pyear >= 1963,]
years <- min(compcrsp$pyear):max(compcrsp$pyear)
min.year <- min(years)
breakpoints <- list()
for (i in 1:length(vars)) {
    var <- vars[i]
    if (annual[i]) {
        breakpoints[[var]] <- data.frame(pyear=years)
        quantiles <- t(apply(breakpoints[[var]], 1, calc.breakpoints, compcrsp, var, c(0.2, 0.8), TRUE))
        colnames(quantiles) <- c('q1.bp', 'q5.bp')
        breakpoints[[var]] <- cbind(breakpoints[[var]], quantiles)
    } else {
        breakpoints[[var]] <- expand.grid(pyear=years, pmonth=1:12)
        quantiles <- t(apply(breakpoints[[var]], 1, calc.breakpoints, compcrsp, var, c(0.2, 0.8), FALSE))
        colnames(quantiles) <- c('q1.bp', 'q5.bp')
        breakpoints[[var]] <- cbind(breakpoints[[var]], quantiles)
    }
}

calc.anom <- function(anom, compcrsp, var) {
    date.indices <- compcrsp$pyear == anom[1] & compcrsp$pmonth == anom[2]
    q1.comb.indices <- which(date.indices & compcrsp[var] <= anom[3])
    q5.comb.indices <- which(date.indices & compcrsp[var] >= anom[4])
    q1.weights <- compcrsp$lagmktcap[q1.comb.indices]
    q5.weights <- compcrsp$lagmktcap[q5.comb.indices]
    q1.returns <- compcrsp$ret[q1.comb.indices]
    q5.returns <- compcrsp$ret[q5.comb.indices]
    q1vwret <- sum(q1.weights/sum(q1.weights, na.rm=TRUE) * q1.returns, na.rm=TRUE)
    q5vwret <- sum(q5.weights/sum(q5.weights, na.rm=TRUE) * q5.returns, na.rm=TRUE)
    q1ewret <- mean(q1.returns, na.rm=TRUE)
    q5ewret <- mean(q5.returns, na.rm=TRUE)
    return(c(q1vwret=q1vwret, q5vwret=q5vwret, q1ewret=q1ewret, q5ewret=q5ewret))
}

# calculate the value-weighted and equally-weighted portfolio returns for each accounting variable
nmonth <- 12 * (max(years) - min(years)) + tail(compcrsp$pmonth, 1)
anom <- list()
for (i in 1:length(vars)) {
    var <- vars[i]
    print(var)
    
    anom.var <- data.frame(pyear=head(rep(years, each=12), nmonth), 
                           pmonth=head(rep(1:12, times=length(years)), nmonth))
    
    if (annual[i]) {
        anom.var <- merge(anom.var, breakpoints[[var]], by=c('pyear'))
    } else {
        anom.var <- merge(anom.var, breakpoints[[var]], by=c('pyear', 'pmonth'))
    }
    
    anom.var <- anom.var[order(anom.var$pyear, anom.var$pmonth),]    
    anom.var <- cbind(anom.var, t(apply(anom.var, 1, calc.anom, compcrsp, var)))
    anom.var$zcvwret <- anom.var$q5vwret - anom.var$q1vwret
    anom.var$zcewret <- anom.var$q5ewret - anom.var$q1ewret
    anom.var$year <- ifelse(anom.var$pmonth > 6, anom.var$pyear + 1, anom.var$pyear)
    anom.var$month <- ifelse(anom.var$pmonth > 6, anom.var$pmonth - 6, anom.var$pmonth + 6)
    
    anom[[var]] <- anom.var
}

# load the Fama-French factor data
factors <- read.csv("data/F-F_Research_Data_Factors.csv", header=TRUE)

# convert the dates (originally in YYYYMM integer form) into years and month
factors$year <- as.numeric(substr(as.character(factors$date), 1, 4))
factors$month <- as.numeric(substr(as.character(factors$date), 5, 6))
factors <- subset(factors, select=c("year", "month", "EXMKT", "SMB", "HML", "RF"))
# calculate pyear and pmonth
factors$pyear <- ifelse(factors$month < 7, factors$year - 1, factors$year)
factors$pmonth <- ifelse(factors$month < 7, factors$month + 6, factors$month - 6)

# initialize arrays for the annual coefficient estimates
capm.zc.vw <- list()
capm.zc.ew <- list()
ff.zc.vw <- list()
ff.zc.ew <- list()
capm <- list()
ff <- list()

summary.stats <- function(x) {
    results <- matrix(0, 4, ncol(x))
    rownames(results) <- c('mean', 'sd', 'se', 't')
    colnames(results) <- colnames(x)
    results['mean',] <- apply(x, 2, mean)
    results['sd',] <- apply(x, 2, sd)
    results['se',] <- results['sd',] / sqrt(length(years))
    results['t',] <- results['mean',] / results['se',]
    return(results)
}

# perform Fama-French regressions
for (var in vars) {
    anom.ff <- merge(anom[[var]], factors, by=c('pyear', 'pmonth'))
    capm.zc.vw[[var]] <- matrix(NaN, length(years), 2)
    capm.zc.ew[[var]] <- matrix(NaN, length(years), 2)
    ff.zc.vw[[var]] <- matrix(NaN, length(years), 4)
    ff.zc.ew[[var]] <- matrix(NaN, length(years), 4)
    colnames(capm.zc.vw[[var]]) <- c('vw.alpha', 'vw.beta')
    colnames(capm.zc.ew[[var]]) <- c('ew.alpha', 'ew.beta')
    colnames(ff.zc.vw[[var]]) <- c('vw.alpha', 'vw.beta', 'vw.smb', 'vw.hml')
    colnames(ff.zc.ew[[var]]) <- c('ew.alpha', 'ew.beta', 'ew.smb', 'ew.hml')
    
    for (t in 1:length(years)) {
        rows <- anom.ff$pyear == years[t]
        y.vw <- anom.ff$zcvwret[rows]
        y.ew <- anom.ff$zcewret[rows]
        x.capm <- cbind(rep(1, length(y.vw)), anom.ff$EXMKT[rows] / 100)
        x.ff <- cbind(rep(1, length(y.vw)),
                   anom.ff$EXMKT[rows] / 100,
                   anom.ff$SMB[rows] / 100,
                   anom.ff$HML[rows] / 100)
        capm.zc.vw[[var]][t,] <- solve(t(x.capm) %*% x.capm, t(x.capm) %*% y.vw)
        capm.zc.ew[[var]][t,] <- solve(t(x.capm) %*% x.capm, t(x.capm) %*% y.ew)
        ff.zc.vw[[var]][t,] <- solve(t(x.ff) %*% x.ff, t(x.ff) %*% y.vw)
        ff.zc.ew[[var]][t,] <- solve(t(x.ff) %*% x.ff, t(x.ff) %*% y.ew)
    }
    
    capm[[var]] <- cbind(summary.stats(capm.zc.vw[[var]]), summary.stats(capm.zc.ew[[var]]))
    ff[[var]] <- cbind(summary.stats(ff.zc.vw[[var]]), summary.stats(ff.zc.ew[[var]]))
}

save(anom, file='anom.RData')
save(capm, file='capm.RData')
save(ff, file='ff.RData')