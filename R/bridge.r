# notes ----
# bridging GOA northern rockfish from ADMB to RTMB
# ben.williams@noaa.gov
# last update
# 2024-09

# load ----
library(RTMB)
library(tidyverse)
library(Matrix)
library(tmbstan)
library(shinystan)
library(here)
library(scico)
theme_set(theme_bw())
# devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
library(adnuts)
# install.packages('StanEstimators', repos = c('https://andrjohns.r-universe.dev', 'https://cloud.r-project.org'))
library(StanEstimators)

source(here::here('r', 'bridge_data.r'))
source(here::here('r', 'bridge_pars.r'))
source(here::here('r', 'bridge_model.r'))
source(here::here('r', 'utils.r'))

# ADMB model output 
catch <- read.csv(here::here('admb', 'processed', 'catch.csv'))
srv <- read.csv(here::here('admb', 'processed', 'survey.csv'))
fac <- read.csv(here::here('admb', 'processed', 'fac.csv'))
sac <- read.csv(here::here('admb', 'processed', 'sac.csv'))
fsc <- read.csv(here::here('admb', 'processed', 'fsc.csv'))
slx <- read.csv(here::here('admb', 'processed', 'selex.csv'))
bio <- read.csv(here::here('admb', 'processed', 'bio_rec_f.csv'))
b40 <- read.csv(here::here('admb', 'processed', 'b35_b40_yld.csv'))
REP <- readLines(here::here('admb', "nr.rep"))
# data ----
# create data list for RTMB
data <- list(
             ages = ages,
             years = years,
             length_bins = 15:45,
             waa = waa,
             maa = maa,
             wt_mature = maa * waa / 2, # females only
             spawn_mo = 5,
             catch_ind = rep(1, length(years)),
             catch_obs = catch_obs,
             catch_wt = c(rep(5, 17), rep(50, 45)),
             srv_obs = srv_obs,
             srv_ind = srv_ind,
             srv_sd = srv_sd,
             srv_wt = 0.25,
             fish_age_obs = fish_age_obs,
             fish_age_ind = fish_age_ind,
             fish_age_iss = fish_age_iss,
             fish_age_wt = 0.5,
             srv_age_obs = srv_age_obs,
             srv_age_ind = srv_age_ind,
             srv_age_iss = srv_age_iss,
             srv_age_wt = 0.5,
             fish_size_obs = fish_size_obs,
             fish_size_ind = fish_size_ind,
             fish_size_iss = fish_size_iss,
             fish_size_wt = 0.5,
             age_error = age_error,
             size_age = size_age,
             wt_fmort_reg = 0.1,
             wt_rec_var = 1,
             mean_M = 0.06,
             cv_M = 0.05,
             mean_q = 1,
             cv_q = 0.45,
             mean_sigmaR = 1.5,
             cv_sigmaR = 0.01,
             yield_ratio = yield_ratio
)

# list of pars from ADMB model
pars = list(log_M = log_M,
        log_a50C = log_a50C,
        deltaC = deltaC,
        log_a50S = log_a50S,
        deltaS = deltaS,
        log_q = log_q,
        log_mean_R = log_mean_R,
        init_log_Rt = init_log_Rt,
        log_Rt = log_Rt,
        log_mean_F = log_mean_F,
        log_Ft =  log_Ft,
        log_F35 = log_F35,
        log_F40 = log_F40,
        log_F50 = log_F50,
        sigmaR = sigmaR)

map = list(log_M = factor(NA),
           log_a50C = factor(NA),
           deltaC = factor(NA),
           log_a50S = factor(NA),
           deltaS = factor(NA),
           log_q = factor(NA),
           log_mean_R = factor(NA),
           init_log_Rt = factor(rep(NA, length(init_log_Rt))),
           log_Rt = factor(rep(NA, length(log_Rt))),
           log_mean_F = factor(NA),
           log_Ft = factor(rep(NA, length(log_Ft))),
           log_F35 = factor(NA),
           log_F40 = factor(NA),
           log_F50 = factor(NA),
           sigmaR = factor(NA))

# build the model 
# comparison run to admb model
obj <- RTMB::MakeADFun(f, 
                        pars, 
                        map = map)  
report <- obj$report(obj$env$last.par.best)
# some comparisons - values match up well
proj_bio(report)
REP[grep('spawn_biom', REP) + 1]
REP[grep('tot_biom', REP) + 1]
REP[grep('ABC', REP)[1] + 1] # F40
REP[grep('ABC', REP)[3] + 1] # ABC catch - 2023
REP[grep('ABC', REP)[4] + 1] # ABC catch - 20234
REP[grep('OFL', REP)[1] + 1] # FOFL
REP[grep('OFL', REP)[3] + 1] # FOFL catch - 2023
REP[grep('OFL', REP)[4] + 1] # FOFL catch - 2024
REP[grep('q_trawl', REP) + 1] # q
report$q
REP[grep('nat_mort', REP) + 1] # q
report$M
REP[grep('B40', REP) + 1]
report$B40
REP[grep('B_zero', REP) + 1] # B0
report$B0

# estimate RTMB model ----
# examine model run with same starting point as ADMB
pars2 = list(log_M = log(0.06),
             log_a50C = log(7.5),
             deltaC = 3.0,
             log_a50S = log(7.3),
             deltaS = 3.8,
             log_q = log(1),
             log_mean_R = 4.3,
             init_log_Rt = rep(0, length(pars$init_log_Rt)),
             log_Rt = rep(0, length(years)),
             log_mean_F = 0,
             log_Ft =  rep(0, length(years)),
             log_F35 = 0,
             log_F40 = 0,
             log_F50 = 0,
             sigmaR = 1.5)

# parameter bounds - same as ADMB
lower = c(-Inf, # log M
          -Inf, #log a50C
          -Inf, # delta C
          -Inf, # log_a50S
          -Inf, # delta S
          -Inf, # logq
          -15, # log mean R
          rep(-10, length(pars$init_log_Rt)), # init rec devs
          rep(-10, length(years)), # rec devs
          -15, # log mean F
          rep(-15, length(years)), # Fdevs
          rep(-4.605,3)) # Fspr
          #0.3) # sigmaR

upper = c(Inf, # log M
          Inf, #log a50C
          Inf, # delta C
          Inf, # log_a50S
          Inf, # delta S
          Inf, # logq
          10, # log mean R
          rep(10,  length(pars$init_log_Rt)), # init rec devs
          rep(10, length(years)), # rec devs
          15, # log mean F
          rep(15, length(years)), # Fdevs
          rep(0,3)) # Fspr
          # 10) # sigmaR

obj1 <- RTMB::MakeADFun(f, 
                        pars2,
                        map = list(sigmaR = factor(NA)))  
# obj1$report(obj1$env$last.par.best)$spawn_bio
fit1 <- nlminb(start = obj1$par,
               objective = obj1$fn,
               gradient = obj1$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
                lower = lower,
                upper = upper)

report1 <- obj1$report(obj1$env$last.par.best)
proj_bio(report1) # these differ from the orig
proj_bio(report) # bridge results
report1$M # M to the floor
report$M # orig M
report1$q # q to the roof!
report$q # orig q

# fixed M ----
pars2$log_M = log(report$M)

# parameter bounds
lower = c(#-Inf, # log M
          -Inf, #log a50C
          -Inf, # delta C
          -Inf, # log_a50S
          -Inf, # delta S
          -Inf, # logq
          -15, # log mean R
          rep(-10, length(pars$init_log_Rt)), # init rec devs
          rep(-10, length(years)), # rec devs
          -15, # log mean F
          rep(-15, length(years)), # Fdevs
          rep(-4.605,3)) # Fspr
          #0.3) # sigmaR

upper = c(#Inf, # log M
          Inf, #log a50C
          Inf, # delta C
          Inf, # log_a50S
          Inf, # delta S
          Inf, # logq
          10, # log mean R
          rep(10,  length(pars$init_log_Rt)), # init rec devs
          rep(10, length(years)), # rec devs
          15, # log mean F
          rep(15, length(years)), # Fdevs
          rep(0,3)) # Fspr
          # 10) # sigmaR


obj2 <- RTMB::MakeADFun(f, 
                        pars2,
                        map = list(log_M = factor(NA),
                                   sigmaR = factor(NA)))  
# obj1$report(obj1$env$last.par.best)$spawn_bio
fit2 <- nlminb(start = obj2$par,
               objective = obj2$fn,
               gradient = obj2$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
               lower = lower,
               upper = upper)

report2 <- obj2$report(obj2$env$last.par.best)
proj_bio(report2) # still different, but closer
proj_bio(report)
report2$q # q is down a bit, so biomass a bit higher
report$q # orig q

# can get SD for all params
sdrep2 <- sdreport(obj2, getJointPrecision = TRUE)

# Profile the likelihood
multi_like(par_names = c("log_a50C", "log_a50S"), 
           par_ranges = list(log_a50C = seq(log(3), log(12), length.out=50),
                             log_a50S = seq(log(3), log(12), length.out=50)),
           obj = obj2, 
           fit =fit2) %>% 
  ggplot(aes(log_a50C, log_a50S)) +
  geom_tile(aes(fill = log_like)) +
  stat_contour(aes(z = log_like)) +
  scico::scale_fill_scico(palette = 'oslo')


logq = seq(log(0.5), log(1.5), length.out=50)
plot(logq,
single_like(par_name = 'log_q',
            par_values = logq,
            obj=obj2,
            fit = fit2), type = 'l')

loga50C = seq(log(5), log(11), length.out=50)
plot(loga50C,
single_like(par_name = 'log_a50C',
            par_values = loga50C,
            obj=obj2,
            fit = fit2), type = 'l')     
loga50S = seq(log(5), log(11), length.out=50)
lines(loga50S,
single_like(par_name = 'log_a50S',
            par_values = loga50S,
            obj=obj2,
            fit = fit2), col=4)                        

Q <- sdrep2$jointPrecision
## if no random effects this will be NULL
if(!is.null(Q)){
  M <- solve(Q)
} else {
  M <- sdrep2$cov.fixed
  Q <- solve(M)
}         

# Cole recommends 5 chains, 1000 iters, 250 warmup
# has substantial divergence
mcmc2 <- sample_sparse_tmb(obj2, iter=1000, warmup=250,
                          init='random', chains=5,
                          cores=3, metric='dense',
                          Qinv=M, Q=Q,
                          globals=list(data = data), skip_optimization=TRUE)

summary(mcmc2$monitor$n_eff)
summary(mcmc2$monitor$Rhat)
plot_marginals(mcmc2, pars=1:5)

# selectity priors ----
# also estimating M
data$mean_a50C = 7.5
data$cv_a50C = 1
data$mean_deltaC = 3.8
data$cv_deltaC = 1
data$mean_q
data$cv_q = 0.45
data$cv_M = 0.05

lower = c(-Inf, # log M
          -Inf, #log a50C
          -Inf, # delta C
          -Inf, # log_a50S
          -Inf, # delta S
          -Inf, # logq
  -15, # log mean R
  rep(-10, length(pars$init_log_Rt)), # init rec devs
  rep(-10, length(years)), # rec devs
  -15, # log mean F
  rep(-15, length(years)), # Fdevs
  rep(-4.605,3)) # Fspr

upper = c(Inf, # log M
          Inf, #log a50C
          Inf, # delta C
          Inf, # log_a50S
          Inf, # delta S
          Inf, # logq
  10, # log mean R
  rep(10,  length(pars$init_log_Rt)), # init rec devs
  rep(10, length(years)), # rec devs
  15, # log mean F
  rep(15, length(years)), # Fdevs
  rep(0,3)) # Fspr


# change model to accept priors
obj3 <- RTMB::MakeADFun(f1, 
                        pars2,
                        map = list(sigmaR = factor(NA)))  
# obj1$report(obj1$env$last.par.best)$spawn_bio
fit3 <- nlminb(start = obj3$par,
               objective = obj3$fn,
               gradient = obj3$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
               lower = lower,
               upper = upper)

sdrep3 <- sdreport(obj3, getJointPrecision = TRUE)
report3 <- obj3$report(obj3$env$last.par.best)
proj_bio(report3) # higher values than orig
proj_bio(report)
report3$M # M up a bit
report$M
report3$q # q down a bit
report$q

Q <- sdrep3$jointPrecision
## if no random effects this will be NULL
if(!is.null(Q)){
  M <- solve(Q)
} else {
  M <- sdrep3$cov.fixed
  Q <- solve(M)
}         

# Cole recommends 5 chains, 1000 iters, 250 warmup
mcmc3 <- sample_sparse_tmb(obj3, iter=1200, warmup=400,
                          init='random', chains=5,
                          cores=5, metric='dense',
                          Qinv=M, Q=Q,
                          globals=list(data = data), 
                          skip_optimization=TRUE,
                          control = list(adapt_delta = 0.99))
summary(mcmc3)
summary(mcmc3$monitor$n_eff)
summary(mcmc3$monitor$Rhat)
post <- extract_samples(mcmc3, as.list=TRUE)
postlist <- coda::mcmc.list(lapply(post, coda::mcmc))
coda::traceplot(postlist)
glimpse(post)
sp <- extract_sampler_params(mcmc3)
glimpse(sp)
plot_marginals(mcmc3, pars=1:6)
adnuts::plot_uncertainties(mcmc3)

# get posterior for derived quantities

get_spawn_bio(obj3, post, iters=50) -> ssb
ssb %>% 
  ggplot(aes(year, spawn_bio, group = sim)) + 
  geom_line(alpha = 0.2)

get_tot_bio(obj3, post, iters=50) -> tb
tb %>% 
  ggplot(aes(year, tot_bio, group = sim)) + 
  geom_line(alpha = 0.2)
  
get_rec(obj3, post, iters=50) -> recs
recs %>% 
  ggplot(aes(year, recruits, group = sim)) + 
  geom_line(alpha = 0.2) +
  coord_cartesian(ylim = c(0, 300))

pairs_admb(mcmc3, pars=1:6, order='slow')
pairs_admb(mcmc3, pars=1:6, order='mismatch')
launch_shinyadmb(mcmc3)

multi_like(par_names = c("log_a50C", "log_a50S"),
           par_ranges = list(log_a50C = seq(log(3), log(12), length.out=50),  # Replace with appropriate range
                             log_a50S = seq(log(3), log(12), length.out=50)), 
           obj = obj3, 
           fit = fit3) %>% 
  ggplot(aes(log_a50C, log_a50S)) +
  geom_tile(aes(fill = log_like)) +
  stat_contour(aes(z = log_like)) +
  scico::scale_fill_scico(palette = 'oslo')

multi_like(par_names = c("deltaS", "log_a50S"),
           par_ranges = list(deltaS = seq(1, 8, length.out=50),  # Replace with appropriate range
                             log_a50S = seq(log(4), log(12), length.out=50)), 
           obj = obj3, 
           fit = fit3) %>% 
  ggplot(aes(deltaS, log_a50S)) +
  geom_tile(aes(fill = log_like)) +
  stat_contour(aes(z = log_like)) +
  scico::scale_fill_scico(palette = 'oslo')

multi_like(par_names = c("log_M", "log_q") , 
           par_ranges = 
             list(log_M = seq(log(0.01), log(0.2), length.out=50),  # Replace with appropriate range
                  log_q = seq(log(.1), log(3), length.out=50)), 
               obj = obj3, 
               fit = fit3) %>% 
  ggplot(aes(log_M, log_q)) +
  geom_tile(aes(fill = log_like)) +
  stat_contour(aes(z = log_like)) +
  scico::scale_fill_scico(palette = 'oslo')

loga50C = seq(log(5), log(11), length.out=50)
plot(loga50C,
single_like(par_name = 'log_a50C',
            par_values = loga50C,
            obj=obj3,
            fit = fit3), type = 'l')     
loga50S = seq(log(5), log(11), length.out=50)
lines(loga50S,
single_like(par_name = 'log_a50S',
            par_values = loga50S,
            obj=obj3,
            fit = fit3), col=4)      

logq = seq(log(0.5), log(1.2), length.out=50)
plot(logq,
single_like(par_name = 'log_q',
            par_values = logq,
            obj=obj3,
            fit = fit3), type='l')   

logm = seq(log(0.005), log(0.2), length.out=50)
plot(logm,
single_like(par_name = 'log_M',
            par_values = logm,
            obj=obj3,
            fit = fit3), type='l')   


# selectivity deltas
par_names <- c("deltaC", "deltaS")  # Replace with your parameter names
par_ranges <- list(
  deltaC = seq(.05, 7, length.out=50),  # Replace with appropriate range
  deltaS = seq(.05, 7, length.out=50)   # Replace with appropriate range
)

# Profile the likelihood
multi_like(par_names, par_ranges, obj2, fit2) %>% 
  ggplot(aes(deltaC, deltaS)) +
  geom_tile(aes(fill = log_like)) +
  stat_contour(aes(z = log_like)) +
  scico::scale_fill_scico(direction = -1, palette = 'oslo')



np = extract_sampler_params(mcmc3) %>%
  pivot_longer(-c(chain, iteration), names_to='Parameter', values_to='Value') %>% 
  select(Iteration=iteration, Parameter, Value, Chain=chain) %>%
  mutate(Parameter=factor(Parameter),
         Iteration=as.integer(Iteration),
         Chain=as.integer(Chain)) 
  


# comparisons ----
# ignore this stuff
## selectivity
as.data.frame(report$slx) %>% 
  rename(fishery = V1, survey = V2) %>% 
  bind_cols(as.data.frame(report3$slx) %>% 
              rename(fishery.1 = V1, survey.1 = V2)) %>% 
  bind_cols(slx) -> slxs 
  
slxs %>% 
  rename(`fishery-rtmb`=fishery, `srv-rtmb`=survey,`fishery-base.1a`=fish,`srv-base.1a`=srv1) %>% 
  pivot_longer(-c(age, maturity)) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_line() +
  scico::scale_color_scico_d(name="", palette = 'roma') +
  theme(legend.position = c(0.8, 0.2))

# ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "slx.png"), width=6.5, height=4.5, units='in')

slxs %>% 
  mutate(fishery = (fishery - fish)/fish,
         survey = (survey - srv1)/srv1) %>% 
  dplyr::select(age, fishery, survey) %>% 
  tidyr::pivot_longer(-age) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  scico::scale_color_scico_d(name="",palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3) +
  theme(legend.position = c(0.8, 0.2))
# ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "slx-diff.png"), width=6.5, height=4.5, units='in')

## biomass
data.frame('ssb-RTMB' = report$spawn_bio,
           'tot-RTMB' = report$tot_bio) %>% 
  bind_cols(data.frame('ssb-RTMB.1' = report3$spawn_bio,
                       'tot-RTMB.1' = report3$tot_bio)) %>% 
  bind_cols(bio) %>% 
  rename(`tot-Base.1a`=tot_biom, `ssb-Base.1a`=sp_biom) %>% 
  pivot_longer(-c(year, F, recruits)) -> Bs
  
Bs %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_line() +
  scico::scale_color_scico_d(name="",palette = 'roma') +
  expand_limits(y=0) +
  theme(legend.position = c(0.8, 0.2)) +
  ylab('Biomass') +
  xlab('Year')
# ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "bio.png"), width=6.5, height=4.5, units='in')


Bs %>% 
  pivot_wider(values_from=value, names_from=name) %>% 
  mutate(tot = (`tot.RTMB` - `tot-Base.1a`) / `tot-Base.1a`,
         ssb = (`ssb.RTMB` - `ssb-Base.1a`) / `ssb-Base.1a`) %>% 
  pivot_longer(c(tot, ssb)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous(labels = scales::percent) +
  scico::scale_color_scico_d(name="", palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)

ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "bio-diff.png"), width=6.5, height=4.5, units='in')

#F

 Bs %>% 
   mutate(id = case_when(grepl('ssb', name) ~ sub('ssb-', "", name),
                         grepl('tot', name) ~ sub('tot-', "", name))) %>% 
   distinct(year, F, id) %>% 
   filter(!(id %in% c('ssb.RTMB', 'ssb.RTMB.1'))) -> Fs
   
 Fs %>% 
   ggplot(aes(year, F, color = id)) +
   geom_line() +
   scico::scale_color_scico_d(name="",palette = 'roma') +
   expand_limits(y=0) +
   theme(legend.position = c(0.8, 0.4)) +
   ylab('F') +
   xlab('Year')

  Fs %>% 
    distinct(year, F, id) %>% 
    pivot_wider(values_from=F, names_from=id) %>% 
    mutate(diff = (tot.RTMB - `Base.1a`) / `Base.1a`) %>% 
  ggplot(aes(year, diff)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)
  
  
# survey biomass
data.frame(RTMB = report$srv_pred) %>% 
  bind_cols(RTMB.1 = report3$srv_pred) %>% 
  bind_cols(srv) %>% 
  rename(Base.1a = pred) -> ds

  ds %>% 
    pivot_longer(c(RTMB, Base.1a, RTMB.1)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point(aes(y=biomass), color = 'darkgray') +
    geom_errorbar(aes(ymin=lci, ymax=uci), color = 'darkgray', width = 0.2) +
  geom_line() +
  scico::scale_color_scico_d(name="",palette='roma', end = 0.6) +
  expand_limits(y=0) +
  ylab('Biomass') +
  xlab('Years') +
  theme(legend.position = c(0.7, 0.2))

  # ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "srvbio.png"), width=6.5, height=4.5, units='in')
  
  
  ds %>% 
    select(year, RTMB, Base.1a) %>% 
    mutate(diff = (RTMB - Base.1a) / Base.1a) %>% 
    ggplot(aes(year, diff)) + 
    geom_point() +
    xlab('Year') +
    scale_y_continuous(labels = scales::percent) +
    # scico::scale_color_scico_d(name="", palette = 'roma') +
    ylab('percent difference') +
    geom_hline(yintercept=0, lty=3)
  
  # ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "srvbio-diff.png"), width=6.5, height=4.5, units='in')
  
## fish age comp
data.frame(age = data$ages, 
           report$fish_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(fish_age_yrs, length(ages)),
         groups = 'rtmb') %>% 
  bind_rows(data.frame(age = data$ages, 
                       report1$fish_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(fish_age_yrs, length(ages)),
                     groups = 'rtmb.1')) %>% 
  bind_rows(fac %>% 
              mutate(groups = 'admb.1a'))-> df5

df5 %>%  
  # mutate(groups = ifelse(groups=='pred', 'admb.1a', groups)) %>% distinct(groups)
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(19,19, 19),guide = guide_none()) +
  # scale_linetype_manual(values = c(1,0,1,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)

### Fishery size compositions

data.frame(age = data$ages, 
           report$fish_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(fish_age_yrs, length(ages)),
         groups = 'RTMB') %>% 
  bind_rows(data.frame(age = data$ages, 
                       report3$fish_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(fish_age_yrs, length(ages)),
                     groups = 'RTMB.1')) %>% 
  bind_rows(fac) -> df5

df5 %>% 
  mutate(groups = ifelse(groups=='pred', 'Base.1a', groups)) %>% 
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(NA,19,NA,NA),guide = guide_none()) +
  scale_linetype_manual(values = c(1,0,1,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)

# ggsave('fsc.png', width=6.5, height=6.5, units='in')

df5 %>% 
  select(age, value, year, groups) %>% 
  pivot_wider(names_from=groups, values_from = value) %>% 
  mutate(diff = (RTMB - pred) / pred,
         Age = factor(age)) %>% 
  pivot_longer(diff) %>% 
  ggplot(aes(year, value, color = Age)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)

# ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "fac-diff.png"), width=6.5, height=6.5, units='in')

# changes weights to 1 --------------
data$srv_wt = 0.25
data$fish_age_wt = 0.5
data$srv_age_wt = 0.5
data$fish_size_wt = 0.5
data$wt_rec_var = 1
data$wt_fmort_reg = 0.1

obj4 <- RTMB::MakeADFun(f1, 
                        pars2,
                        map = list(sigmaR = factor(NA)))  
# obj1$report(obj1$env$last.par.best)$spawn_bio
fit4 <- nlminb(start = obj4$par,
               objective = obj4$fn,
               gradient = obj4$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
               lower = lower,
               upper = upper)

sdrep4 <- sdreport(obj4, getJointPrecision = TRUE)
report4 <- obj4$report(obj4$env$last.par.best)
proj_bio(report4) # higher values than orig
proj_bio(report3)
report4$M # M up a bit
report$M
report4$q # q down a bit
report$q

report4$nll
report3$nll


Q <- sdrep3$jointPrecision
## if no random effects this will be NULL
if(!is.null(Q)){
  M <- solve(Q)
} else {
  M <- sdrep3$cov.fixed
  Q <- solve(M)
}         

# Cole recommends 5 chains, 1000 iters, 250 warmup
mcmc3 <- sample_sparse_tmb(obj3, iter=1200, warmup=400,
                          init='random', chains=5,
                          cores=5, metric='dense',
                          Qinv=M, Q=Q,
                          globals=list(data = data), 
                          skip_optimization=TRUE,
                          control = list(adapt_delta = 0.99))
summary(mcmc3)
summary(mcmc3$monitor$n_eff)
summary(mcmc3$monitor$Rhat)
post <- extract_samples(mcmc3, as.list=TRUE)
postlist <- coda::mcmc.list(lapply(post, coda::mcmc))
coda::traceplot(postlist)
glimpse(post)
sp <- extract_sampler_params(mcmc3)
glimpse(sp)
plot_marginals(mcmc3, pars=1:6)
adnuts::plot_uncertainties(mcmc3)

# get posterior for derived quantities

get_spawn_bio(obj3, post, iters=50) -> ssb
ssb %>% 
  ggplot(aes(year, spawn_bio, group = sim)) + 
  geom_line(alpha = 0.2)

get_tot_bio(obj3, post, iters=50) -> tb
tb %>% 
  ggplot(aes(year, tot_bio, group = sim)) + 
  geom_line(alpha = 0.2)
  
get_rec(obj3, post, iters=50) -> recs
recs %>% 
  ggplot(aes(year, recruits, group = sim)) + 
  geom_line(alpha = 0.2) +
  coord_cartesian(ylim = c(0, 300))

pairs_admb(mcmc3, pars=1:6, order='slow')
pairs_admb(mcmc3, pars=1:6, order='mismatch')
launch_shinyadmb(mcmc3)