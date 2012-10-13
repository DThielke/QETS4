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
returns <- crsp.clean[, c('PERMNO', 'pyear', 'pmonth', 'mktcap', 'RET', 'EXCHCD')]
names(returns) <- tolower(names(returns))

# merge the accounting and returns data
compcrsp <- merge(comp, returns, by=c('permno', 'pyear'), all.y=TRUE)
compcrsp <- compcrsp[compcrsp$pyear %in% comp$pyear,]
compcrsp <- compcrsp[order(compcrsp$permno, compcrsp$pyear, compcrsp$pmonth),]

# accounting variables
vars <- c('btm', 'roa', 'agr', 'nsi', 'acc')

# calculate breakpoints for BTM, ROA, asset growth, net stock issues, and accruals
years <- min(compcrsp$pyear):max(compcrsp$pyear)
probabilities <- c(0.2, 0.4, 0.6, 0.8)
breakpoints <- list()
for (var in vars) {
    breakpoints[[var]] <- matrix(0, length(years), 4)
    for (t in 1:length(years)) {
        indices <- compcrsp$exchcd == 1 & compcrsp$pyear == years[t]
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
    pyear <- min(years)
    pmonth <- 1
    for (t in 1:nmonth) {
        date.indices <- compcrsp$pyear == pyear & compcrsp$pmonth == pmonth
        q1.indices <- compcrsp[var] <= breakpoints[[var]][floor(t/12) + 1, 1]
        q5.indices <- compcrsp[var] > breakpoints[[var]][floor(t/12) + 1, 4]
        q1.weights <- compcrsp$mktcap[which(date.indices & q1.indices)]
        q5.weights <- compcrsp$mktcap[which(date.indices & q5.indices)]
        q1.returns <- compcrsp$ret[which(date.indices & q1.indices)]
        q5.returns <- compcrsp$ret[which(date.indices & q5.indices)]
        anom[[var]][t,'q1vwret'] <- sum(q1.weights/sum(q1.weights) * q1.returns)
        anom[[var]][t,'q5vwret'] <- sum(q5.weights/sum(q5.weights) * q5.returns)
        anom[[var]][t,'q1ewret'] <- mean(q1.returns)
        anom[[var]][t,'q5ewret'] <- mean(q5.returns)
        pmonth <- pmonth + 1
        if (pmonth > 12) {
            pmonth <- 1
            pyear <- pyear + 1
        }
    }
    anom[[var]][is.nan(anom[[var]][,'q1ewret']), 'q1ewret'] <- 0
    anom[[var]][is.nan(anom[[var]][,'q5ewret']), 'q5ewret'] <- 0
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

# initialize arrays for the annual coefficient estimates
ff.zc.vw <- list()
ff.zc.ew <- list()
ff.zc.vw.summary <- list()
ff.zc.ew.summary <- list()

# perform Fama-French regressions
for (var in vars) {
    anom.ff <- merge(anom[[var]], factors, by=c('year', 'month'))
    ff.zc.vw[[var]] <- matrix(NaN, length(years), 4)
    ff.zc.ew[[var]] <- matrix(NaN, length(years), 4)
    colnames(ff.zc.vw[[var]]) <- c('alpha', 'beta', 'gamma', 'delta')
    colnames(ff.zc.ew[[var]]) <- c('alpha', 'beta', 'gamma', 'delta')
    for (t in 1:length(years)) {
        rows <- anom.ff[,'year'] == years[t]
        y.vw <- anom.ff[rows, 'zcvwret']
        y.ew <- anom.ff[rows, 'zcewret']
        x <- cbind(rep(1, length(y.vw)),
                   anom.ff$EXMKT[rows] / 100,
                   anom.ff$SMB[rows] / 100,
                   anom.ff$HML[rows] / 100)
        ff.zc.vw[[var]][t,] <- solve(t(x) %*% x, t(x) %*% y.vw)
        ff.zc.ew[[var]][t,] <- solve(t(x) %*% x, t(x) %*% y.ew)
    }
    
    # value weighted summary statistics
    ff.zc.vw.summary[[var]] <- matrix(0, 4, 4, dimnames=list(c('mean', 'sd', 'se', 't'), c('alpha', 'beta', 'gamma', 'delta')))
    ff.zc.vw.summary[[var]]['mean',] <- apply(ff.zc.vw[[var]], 2, mean)
    ff.zc.vw.summary[[var]]['sd',] <- apply(ff.zc.vw[[var]], 2, sd)
    ff.zc.vw.summary[[var]]['se',] <- sd / sqrt(length(years))
    ff.zc.vw.summary[[var]]['t',] <- sd / se
    
    # equally weighted summary statistics
    ff.zc.ew.summary[[var]] <- matrix(0, 4, 4, dimnames=list(c('mean', 'sd', 'se', 't'), c('alpha', 'beta', 'gamma', 'delta')))
    ff.zc.ew.summary[[var]]['mean',] <- apply(ff.zc.ew[[var]], 2, mean)
    ff.zc.ew.summary[[var]]['sd',] <- apply(ff.zc.ew[[var]], 2, sd)
    ff.zc.ew.summary[[var]]['se',] <- sd / sqrt(length(years))
    ff.zc.ew.summary[[var]]['t',] <- sd / se
}

save(anom, file='anom.RData')
save(ff.zc.vw.summary, file='ff.zc.vw.summary.RData')
save(ff.zc.ew.summary, file='ff.zc.ew.summary.RData')