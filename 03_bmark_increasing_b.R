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

x200 = simulate_data(lmark = 2, n = 200, p = 10)

nrepl = 10
n.cores = 1
time_exp3 = benchmark(
  'n = 200, p = 10, b = 50, nc = 1' = run_prc(x200, 50, n.cores), 
  'n = 200, p = 10, b = 100, nc = 1' = run_prc(x200, 100, n.cores), 
  'n = 200, p = 10, b = 200, nc = 1' = run_prc(x200, 200, n.cores), 
  'n = 200, p = 10, b = 300, nc = 1' = run_prc(x200, 300, n.cores), 
  'n = 200, p = 10, b = 400, nc = 1' = run_prc(x200, 400, n.cores), 
  'n = 200, p = 10, b = 500, nc = 1' = run_prc(x200, 500, n.cores), 
  replications = nrepl
)

save(time_exp3, nrepl, file = 'results/03_results.Rdata')
