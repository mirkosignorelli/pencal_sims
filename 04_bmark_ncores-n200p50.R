library(pencal)
library(rbenchmark)

set.seed(19931101)
n = 200
p = 50
B = 50
nrepl = 10

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

df = simulate_data(lmark = 2, n = n, p = p)

# prerun steps 1 and 2 as input for the timing of steps 2 and 3
y.names = paste('marker', 1:df$p, sep = '')
s1 = fit_lmms(y.names = y.names, 
                 fixefs = ~ age, ranefs = ~ age | id, 
                 long.data = df$simdata$long.data, 
                 surv.data = df$simdata$surv.data,
                 t.from.base = t.from.base,
                 n.boots = B, n.cores = 2, 
                 verbose = F)
s2 = summarize_lmms(object = s1, n.cores = 2, verbose = F)

# functions to replicate steps 1, 2 and 3
step1 = function(df, n.cores) {
  fit_lmms(y.names = y.names, 
                fixefs = ~ age, ranefs = ~ age | id, 
                long.data = df$simdata$long.data, 
                surv.data = df$simdata$surv.data,
                t.from.base = t.from.base,
                n.boots = B, n.cores = n.cores, 
                verbose = T)
}

step2 = function(df, n.cores) {
  summarize_lmms(object = s1, n.cores = n.cores, verbose = T)
}

step3 = function(df, n.cores, pen) {
  fit_prclmm(object = s2, surv.data = df$simdata$surv.data,
             baseline.covs = ~ baseline.age,
             penalty = pen, n.cores = n.cores, verbose = T)
}

# benchmarking step1
t_step1 = benchmark(
  'nc = 1' = step1(df, 1), 
  'nc = 2' = step1(df, 2), 
  'nc = 3' = step1(df, 3), 
  'nc = 4' = step1(df, 4),
  'nc = 8' = step1(df, 8), 
  'nc = 16' = step1(df, 16),
  replications = nrepl
)

# benchmarking step2
t_step2 = benchmark(
  'nc = 1' = step2(df, 1), 
  'nc = 2' = step2(df, 2), 
  'nc = 3' = step2(df, 3), 
  'nc = 4' = step2(df, 4),
  'nc = 8' = step2(df, 8), 
  'nc = 16' = step2(df, 16),
  replications = nrepl
)

# benchmarking step3
# ridge penalty
t_step3_ridge = benchmark(
  'nc = 1, ridge' = step3(df, 1, 'ridge'), 
  'nc = 2, ridge' = step3(df, 2, 'ridge'), 
  'nc = 3, ridge' = step3(df, 3, 'ridge'), 
  'nc = 4, ridge' = step3(df, 4, 'ridge'), 
  'nc = 8, ridge' = step3(df, 8, 'ridge'), 
  'nc = 16, ridge' = step3(df, 16, 'ridge'),
  replications = nrepl
)

# elasticnet penalty
t_step3_elasticnet = benchmark(
  'nc = 1, elasticnet' = step3(df, 1, 'elasticnet'), 
  'nc = 2, elasticnet' = step3(df, 2, 'elasticnet'), 
  'nc = 3, elasticnet' = step3(df, 3, 'elasticnet'), 
  'nc = 4, elasticnet' = step3(df, 4, 'elasticnet'), 
  'nc = 8, elasticnet' = step3(df, 8, 'elasticnet'), 
  'nc = 16, elasticnet' = step3(df, 16, 'elasticnet'),
  replications = nrepl
)

save(n, p, B, nrepl,
     t_step1, t_step2,
     t_step3_ridge, t_step3_elasticnet,
     file = 'results/04_results.Rdata')
