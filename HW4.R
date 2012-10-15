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
returns <- crsp.clean[, c('PERMNO', 'pyear', 'pmonth', 'lagmktcap', 'RET', 'EXCHCD')]
names(returns) <- tolower(names(returns))

# merge the accounting and returns data
compcrsp <- merge(comp, returns, by=c('permno', 'pyear'))
compcrsp <- compcrsp[compcrsp$pyear %in% comp$pyear,]
compcrsp <- compcrsp[order(compcrsp$permno, compcrsp$pyear, compcrsp$pmonth),]

# accounting variables
vars <- c('btm', 'roa', 'agr', 'nsi', 'acc')

# take log of net stock issues
compcrsp$nsi <- log(compcrsp$nsi)

# calculate breakpoints for BTM, ROA, asset growth, net stock issues, and accruals
compcrsp <- compcrsp[compcrsp$pyear >= 1963,]
years <- min(compcrsp$pyear):max(compcrsp$pyear)
probabilities <- c(0.2, 0.8)
breakpoints <- list()
for (var in vars) {
    breakpoints[[var]] <- matrix(0, length(years), 4)
    for (t in 1:length(years)) {
        indices <- compcrsp$exchcd == 1 & compcrsp$pyear == years[t] & compcrsp$pmonth == 1
        breakpoints[[var]][t,] <- quantile(compcrsp[indices, var], probabilities, na.rm=TRUE)
    }
}

# calculate the value-weighted and equally-weighted portfolio returns for each accounting variable
nmonth <- 12 * (max(years) - min(years)) + tail(compcrsp$pmonth, 1)
anom <- list()
for (var in vars) {
    anom[[var]] <- matrix(0, nmonth, 10)
    anom[[var]][,1] <- head(rep(years, each=12), nmonth) # pyear
    anom[[var]][,2] <- head(rep(1:12, times=length(years)), nmonth) # pmonth
    colnames(anom[[var]]) <- c('pyear', 'pmonth', 'year', 'month', 'q1vwret', 'q5vwret', 'q1ewret', 'q5ewret', 'zcvwret', 'zcewret')
    min.year <- min(years)
    pyear <- min.year
    pmonth <- 1
    for (t in 1:nmonth) {
        date.indices <- compcrsp$pyear == pyear & compcrsp$pmonth == pmonth
        q1.indices <- compcrsp[var] <= breakpoints[[var]][pyear - min.year + 1, 1]
        q5.indices <- compcrsp[var] >= breakpoints[[var]][pyear - min.year + 1, 2]
        q1.comb.indices <- which(date.indices & q1.indices)
        q5.comb.indices <- which(date.indices & q5.indices)
        q1.weights <- compcrsp$lagmktcap[q1.comb.indices]
        q5.weights <- compcrsp$lagmktcap[q5.comb.indices]
        q1.returns <- compcrsp$ret[q1.comb.indices]
        q5.returns <- compcrsp$ret[q5.comb.indices]
        anom[[var]][t,'q1vwret'] <- sum(q1.weights/sum(q1.weights, na.rm=TRUE) * q1.returns, na.rm=TRUE)
        anom[[var]][t,'q5vwret'] <- sum(q5.weights/sum(q5.weights, na.rm=TRUE) * q5.returns, na.rm=TRUE)
        anom[[var]][t,'q1ewret'] <- mean(q1.returns, na.rm=TRUE)
        anom[[var]][t,'q5ewret'] <- mean(q5.returns, na.rm=TRUE)
        pmonth <- pmonth + 1
        if (pmonth > 12) {
            pmonth <- 1
            pyear <- pyear + 1
        }
    }
    anom[[var]][,'zcvwret'] <- anom[[var]][,'q5vwret'] - anom[[var]][,'q1vwret']
    anom[[var]][,'zcewret'] <- anom[[var]][,'q5ewret'] - anom[[var]][,'q1ewret']
    anom[[var]][,'year'] <- ifelse(anom[[var]][,'pmonth'] > 6, anom[[var]][,'pyear'] + 1, anom[[var]][,'pyear'])
    anom[[var]][,'month'] <- ifelse(anom[[var]][,'pmonth'] > 6, anom[[var]][,'pmonth'] - 6, anom[[var]][,'pmonth'] + 6)
}

# load the Fama-French factor data
factors <- read.csv("data/F-F_Research_Data_Factors.csv", header=TRUE)

# convert the dates (originally in YYYYMM integer form) into years and month
factors$year <- as.numeric(substr(as.character(factors$date), 1, 4))
factors$month <- as.numeric(substr(as.character(factors$date), 5, 6))
factors <- subset(factors, select=c("year", "month", "EXMKT", "SMB", "HML", "RF"))
factors$pyear <- ifelse(factors$month < 7, factors$year - 1, factors$year)
factors$pmonth <- ifelse(factors$month < 7, factors$month + 6, factors$month - 6)

# initialize arrays for the annual coefficient estimates
capm.zc.vw <- list()
capm.zc.ew <- list()
capm.zc.vw.summary <- list()
capm.zc.ew.summary <- list()
ff.zc.vw <- list()
ff.zc.ew <- list()
ff.zc.vw.summary <- list()
ff.zc.ew.summary <- list()

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
    colnames(capm.zc.vw[[var]]) <- c('alpha', 'beta')
    colnames(capm.zc.ew[[var]]) <- c('alpha', 'beta')
    colnames(ff.zc.vw[[var]]) <- c('alpha', 'beta', 'gamma', 'delta')
    colnames(ff.zc.ew[[var]]) <- c('alpha', 'beta', 'gamma', 'delta')
    
    for (t in 1:length(years)) {
        rows <- anom.ff[,'pyear'] == years[t]
        y.vw <- anom.ff[rows, 'zcvwret']
        y.ew <- anom.ff[rows, 'zcewret']
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
    
    capm.zc.vw.summary[[var]] <- summary.stats(capm.zc.vw[[var]])
    capm.zc.ew.summary[[var]] <- summary.stats(capm.zc.ew[[var]])
    ff.zc.vw.summary[[var]] <- summary.stats(ff.zc.vw[[var]])
    ff.zc.ew.summary[[var]] <- summary.stats(ff.zc.ew[[var]])
}

save(anom, file='anom.RData')
save(ff.zc.vw.summary, file='ff.zc.vw.summary.RData')
save(ff.zc.ew.summary, file='ff.zc.ew.summary.RData')