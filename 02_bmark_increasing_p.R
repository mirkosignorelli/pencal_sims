on.alice = T
install = F

if (on.alice & install) {
  install.packages('/data1/signorellim1/tarballs/pencal_2.0.2.tar.gz', repos = NULL)
  #install.packages('pencal', repos = "https://cloud.r-project.org")
  install.packages('rbenchmark', repos = "https://cloud.r-project.org")
}

if (on.alice) setwd('/data1/signorellim1/06_pencal_benchmark')
if (!on.alice) {
  library(msigno)
  set.cf()
}

library(pencal)
library(rbenchmark)

set.seed(19931101)

simulate_data = function(lmark, n, p) {
  df = simulate_prclmm_data(n = n, p = p, p.relev = ceiling(0.5*p), 
                            t.values = seq(0, lmark, by = 0.5), landmark = lmark,
                            seed = 1, lambda = 0.2, nu = 2,
                            cens.range = c(lmark, 10), 
                            base.age.range = c(3, 5), tau.age = 0.2)
  df$surv.data$time = df$surv.data$time + lmark
  out = list('simdata' = df, 'lmark' = lmark, 'n' = n, 'p' = p)
  return(out)
}

run_prc = function(df, n.boots, n.cores, verb = F) {
  y.names = paste('marker', 1:df$p, sep = '')
  step1 = fit_lmms(y.names = y.names, 
                   fixefs = ~ age, ranefs = ~ age | id, 
                   long.data = df$simdata$long.data, 
                   surv.data = df$simdata$surv.data,
                   t.from.base = t.from.base,
                   n.boots = n.boots, n.cores = n.cores, 
                   verbose = verb)
  step2 = summarize_lmms(object = step1, n.cores = n.cores, verbose = verb)
  step3 = fit_prclmm(object = step2, surv.data = df$simdata$surv.data,
                     baseline.covs = ~ baseline.age,
                     penalty = 'ridge', n.cores = n.cores, verbose = verb)
}

x5 = simulate_data(lmark = 2, n = 200, p = 5)
x10 = simulate_data(lmark = 2, n = 200, p = 10)
x20 = simulate_data(lmark = 2, n = 200, p = 20)
x30 = simulate_data(lmark = 2, n = 200, p = 30)
x40 = simulate_data(lmark = 2, n = 200, p = 40)
x50 = simulate_data(lmark = 2, n = 200, p = 50)

nrepl = 10
n.cores = 1
time_exp2 = benchmark(
  'n = 200, p = 5, b = 0, nc = 1' = run_prc(x5, 0, n.cores),
  'n = 200, p = 10, b = 0, nc = 1' = run_prc(x10, 0, n.cores),
  'n = 200, p = 20, b = 0, nc = 1' = run_prc(x20, 0, n.cores),
  'n = 200, p = 30, b = 0, nc = 1' = run_prc(x30, 0, n.cores),
  'n = 200, p = 40, b = 0, nc = 1' = run_prc(x40, 0, n.cores),
  'n = 200, p = 50, b = 0, nc = 1' = run_prc(x50, 0, n.cores),
  replications = nrepl
)

save(time_exp2, nrepl, file = 'results/02_results.Rdata')
