# Study the dependence of calibration metrics on the bin number
#
varBinSize = function(E, uE, aux) {

  nBins = c(1,2,5,seq(10,160,by=10))
  sel = length(E) / nBins > 30
  nBins = nBins[sel]

  etab = ztab = c()
  for(i in seq_along(nBins)) {
    nBin = nBins[i]
    res = ErrViewLib::plotRelDiag(
      uE, E, ordX = uE, aux = aux,
      nBin = nBin,
      equiPop = TRUE,
      method = 'cho',
      plot = FALSE)
    etab[i] = res$ENCE
    res = ErrViewLib::plotLZV(
      uE, E / uE, aux = aux,
      nBin = nBin,
      equiPop = TRUE,
      method = 'cho',
      score = TRUE,
      plot = FALSE)
    ztab[i] = res$ZVE
  }
  return(
    list(
      nBins = nBins,
      ENCE  = etab,
      ZVE   = ztab
    )
  )
}

figDir = '../Figs'
tabDir = '../Tabs'
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')

# Analytical analysis of ENCE vs bin-size ####

set.seed(123)
M = 5000
u = rep(1,M)
nBins = 2^(0:7)
X = nBins^0.5

# Full Monte Carlo simulation
nTry = 50
ztab = matrix(NA, nrow=nTry, ncol = length(nBins))
for(j in 1:nTry) {
  E = rnorm(M,0,1)
  for(i in seq_along(nBins)) {
    nBin = nBins[i]
    res = ErrViewLib::plotRelDiag(
      u, E,
      nBin = nBin,
      equiPop = TRUE,
      method = 'cho',
      plot = FALSE)
    ztab[j,i] = res$ENCE
  }
}

# Plot and semi-analytical simulation
png(
  file = file.path(figDir, paste0('Fig_01.png')),
  width  = gPars$reso,
  height = gPars$reso
)
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))
par(
  mfrow = c(1, 1),
  mar = mar,
  mgp = mgp,
  pty = 'm',
  tcl = tcl,
  cex = cex,
  lwd = lwd,
  cex.main = 1
)
bias =  c(1, 1.05, 1.1, 1.15, 1.20, 1.25)
for (j in seq_along(bias)) {
  anal = c()
  for(i in seq_along(nBins)) {
    nBin = nBins[i]
    k = round(M/nBin)
    anal[i] = 0
    for (l in 1:nBin) {
      sam = adehabitatLT::rchi(1e5, k) / sqrt(k)
      anal[i] = anal[i] + mean(abs(sam-bias[j])/bias[j])
    }
    anal[i] = anal[i] / nBin
  }
  if(j==1) {
    plot(X, anal, pch = 19, type = 'b', col = 1,
         xlim = c(0, max(X)+2), xaxs = 'i',
         xlab = expression(N^{1/2}),
         ylim = c(0.,0.23),
         ylab = 'ENCE'
    )
    grid()
    matpoints(X, t(ztab),
              pch = 16, col = gPars$cols_tr[5])
    points(X, anal, type = 'b', pch = 19, col = 1)

  } else {
    points(X, anal, type = 'b', pch = 19, col = 1)
  }
  text(X[length(X)], anal[length(X)], round(bias[j],2),
       pos = 4, col = 2, cex=0.75)
  abline(h=abs(1-bias[j])/bias[j], lty = 2, col = 4)
}
box()
dev.off()

regTab = list()
# BUS2022 ####
set = 'QM9'
D = read.table(
  '../Data/BUS2022/qm9_U0_test_Orig.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

S  = D[, "formula"]
R  = D[, "U0"]
C  = D[, "prediction"]
uC = D[, "uncertainty_total"] ^ 0.5
E  = R - C
M = length(E)
aux = 1:M
uE = uC

resENCE = varBinSize(E, uE, aux)

png(
  file = file.path(figDir, paste0('Fig_02.png')),
  width  = 2*gPars$reso,
  height = gPars$reso
)
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))
par(
  mfrow = c(1, 2),
  mar = mar,
  mgp = mgp,
  pty = 'm',
  tcl = tcl,
  cex = cex,
  lwd = lwd,
  cex.main = 1
)
X = sqrt(resENCE$nBins)
Y = resENCE$ENCE
plot(
  X, Y,
  type ='b', lty = 3,
  pch = 16, col = cols[1],
  xlim = c(0, max(X) + 0.2), xaxs = 'i', xlab = expression(N^{1/2}),
  ylim = c(0, max(Y)), ylab = 'ENCE',
  main = set
)
grid()
abline(h = 0, lty = 2, col = 1)
sel = X > 4
Xs = X[sel]
Ys = Y[sel]
reg = lm(Ys ~ Xs)
regTab$ENCE[[set]] = summary(reg)$coefficients
int = summary(reg)$coefficients[1,1]
Uint = 2*summary(reg)$coefficients[1,2]
abline(reg = reg, lty = 1, lwd = 2, col = cols[1])
segments(
  0.05,int-Uint,
  0.05,int+Uint,
  lwd = 8,
  lty = 1,
  lend = 2,
  col = cols[1])
box()

Y = resENCE$ZVE
plot(
  X, Y,
  type ='b', lty = 3,
  pch = 16, col = cols[1],
  xlim = c(0,max(X) + 0.2), xaxs = 'i', xlab = expression(N^{1/2}),
  ylim = c(1,max(Y)), ylab = 'ZVE'
)
grid()
abline(h = 1, lty = 2, col = 1)
reg = lm(Y ~ X)
regTab$ZVE[[set]] = summary(reg)$coefficients
int = summary(reg)$coefficients[1,1]
Uint = 2*summary(reg)$coefficients[1,2]
abline(reg = reg, lty = 1, lwd = 2, col = cols[1])
segments(
  0.05,int-Uint,
  0.05,int+Uint,
  lwd = 8,
  lty = 1,
  lend = 2,
  col = cols[1])
box()

dev.off()

# PAL2022 ####

sets = c('Diffusion', 'Perovskite')
method = 'RF'
limFit = c(4, 2); names(limFit) = sets
iFig = 2
for (iset in seq_along(sets)) {
  set = sets[iset]
  fRoot = paste0(set, '_', method)

  # Get data
  case = 'Test_uncal'
  fName = file.path('..', 'Data', 'PAL2022', paste0(fRoot, '_', case, '.csv'))
  if (!file.exists(fName))
    next
  print(fRoot)
  Muncal = read.csv(fName, header = TRUE)
  sel = Muncal$uE > 1e-6 * sd(Muncal$E)
  if (sum(!sel) > 0) {
    print(c(case, ' non-pos unc:', sum(!sel)))
    Muncal = Muncal[sel,]
  }
  case = 'Test_cal'
  fName = file.path('..', 'Data', 'PAL2022', paste0(fRoot, '_', case, '.csv'))
  if (!file.exists(fName))
    next
  Mcal = read.csv(fName, header = TRUE)
  sel = Mcal$uE > 1e-6 * sd(Mcal$E)
  if (sum(!sel) > 0) {
    print(c(case, ' non-pos unc:', sum(!sel)))
    Mcal = Mcal[sel,]
  }

  # Generare stats
  resUncal = varBinSize(Muncal$E, Muncal$uE, aux = 1:nrow(Muncal))
  resCal   = varBinSize(Mcal$E, Mcal$uE, aux = 1:nrow(Mcal))

  # Plot
  png(
    file = file.path(figDir, paste0('Fig_0',iFig+iset,'.png')),
    width  = 2*gPars$reso,
    height = gPars$reso
  )
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))
  par(
    mfrow = c(1, 2),
    mar = mar,
    mgp = mgp,
    pty = 'm',
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    cex.main = 1
  )
  X = sqrt(resUncal$nBins)
  Y = resUncal$ENCE
  plot(
    X, Y,
    type ='b', lty = 3,
    pch = 16, col = cols[1],
    xlim = c(0,max(X)+0.2), xaxs = 'i', xlab = expression(N^{1/2}),
    ylim = c(0,max(Y)), ylab = 'ENCE',
    main = set
  )
  grid()
  abline(h = 0, lty = 2, col = 1)
  X = sqrt(resCal$nBins)
  Y = resCal$ENCE
  points(X, Y, type = 'b', lty = 3, pch = 17, col = cols[2])
  sel = X > limFit[set]
  Xs = X[sel]
  Ys = Y[sel]
  reg = lm(Ys ~ Xs)
  regTab$ENCE[[set]] = summary(reg)$coefficients
  int = summary(reg)$coefficients[1,1]
  Uint = 2*summary(reg)$coefficients[1,2]
  abline(reg = reg, lty = 1, lwd = 2, col = cols[2])
  segments(
    0.05,int-Uint,
    0.05,int+Uint,
    lwd = 8,
    lty = 1,
    lend = 2,
    col = cols[2])
  legend(
    'bottomright', bty = 'n',
    legend = c('Uncalibrated','Calibrated'),
    pch = 16:17, col = cols[1:2], lty = 3
  )
  box()

  X = sqrt(resUncal$nBins)
  Y = resUncal$ZVE
  plot(
    X, Y,
    type ='b', lty = 3,
    pch = 16, col = cols[1],
    xlim = c(0,max(X)+0.2), xaxs = 'i', xlab = expression(N^{1/2}),
    ylim = c(1,max(Y)), ylab = 'ZVE'
  )
  grid()
  abline(h = 1, lty = 2, col = 1)
  X = sqrt(resCal$nBins)
  Y = resCal$ZVE
  points(X, Y, type = 'b', lty = 3, pch = 17, col = cols[2])
  sel = X > limFit[set]
  Xs = X[sel]
  Ys = Y[sel]
  reg = lm(Ys ~ Xs)
  regTab$ZVE[[set]] = summary(reg)$coefficients
  int = summary(reg)$coefficients[1,1]
  Uint = 2*summary(reg)$coefficients[1,2]
  abline(reg = reg, lty = 1, lwd = 2, col = cols[2])
  segments(
    0.05,int-Uint,
    0.05,int+Uint,
    lwd = 8,
    lty = 1,
    lend = 2,
    col = cols[2])
  box()
  dev.off()
}


# Regression coefficients table ####
tabENCE = tabZVE = data.frame( Set = NA, Intercept = NA, Slope = NA)
for(set in names(regTab$ENCE)) {
  sum = regTab$ENCE[[set]]
  int = ErrViewLib::prettyUnc(sum[1,1], sum[1,2], numDig = 1)
  slo = ErrViewLib::prettyUnc(sum[2,1], sum[2,2], numDig = 1)
  tabENCE = rbind(
    tabENCE,
    data.frame( Set = set, Intercept = int, Slope = slo)
  )
  sum = regTab$ZVE[[set]]
  int = ErrViewLib::prettyUnc(sum[1,1], sum[1,2], numDig = 1)
  slo = ErrViewLib::prettyUnc(sum[2,1], sum[2,2], numDig = 1)
  tabZVE = rbind(
    tabZVE,
    data.frame( Set = set, Intercept = int, Slope = slo)
  )
}
tabENCE = tabENCE[-1,]
tabZVE  = tabZVE[-1,]
tab = rbind(tabENCE, tabZVE)
sink(file =  file.path(tabDir,'tabFits.tex'))
print(knitr::kable(tab, 'latex'))
sink()
