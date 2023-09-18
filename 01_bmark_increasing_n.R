library(pencal)
library(rbenchmark)

nrepl = 10
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

x100 = simulate_data(lmark = 2, n = 100, p = 10)
x200 = simulate_data(lmark = 2, n = 200, p = 10)
x400 = simulate_data(lmark = 2, n = 400, p = 10)
x600 = simulate_data(lmark = 2, n = 600, p = 10)
x800 = simulate_data(lmark = 2, n = 800, p = 10)
x1000 = simulate_data(lmark = 2, n = 1000, p = 10)

n.cores = 1
time_exp1 = benchmark(
  # no bootstrap, single core
  'n = 100, p = 10, b = 0, nc = 1' = run_prc(x100, 0, n.cores),
  'n = 200, p = 10, b = 0, nc = 1' = run_prc(x200, 0, n.cores), 
  'n = 400, p = 10, b = 0, nc = 1' = run_prc(x400, 0, n.cores),
  'n = 600, p = 10, b = 0, nc = 1' = run_prc(x600, 0, n.cores),
  'n = 800, p = 10, b = 0, nc = 1' = run_prc(x800, 0, n.cores),
  'n = 1000, p = 10, b = 0, nc = 1' = run_prc(x1000, 0, n.cores),
  replications = nrepl
)

save(time_exp1, nrepl, file = 'results/01_results.Rdata')
