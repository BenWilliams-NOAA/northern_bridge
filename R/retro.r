# retrospective functions
peels <- function(model, data, pars, map, peel) {
    p = pars2
    n = length(data$years)
    yr = n-peel
    # update data 
    data$years = head(data$years, -peel) # sets the indexing for T
    data$catch_obs = head(data$catch_obs, -peel)# sets the indexing for T
    data$srv_ind[yr:n] <- 0              # survey biomass index
    data$fish_age_ind[(yr-2):n] <- 0     # fishery age comp index
    data$srv_age_ind[(yr-1):n] <- 0      # survey age comp index
    data$fish_size_ind[(yr-1):n] <- 0    # fishery size comp index

    # update pars
    p$log_Rt = head(p$log_Rt, -peel)
    p$log_Ft = head(p$log_Ft, -peel)
    upr = lapply(p, function(x) x + Inf)
    lwr = lapply(p, function(x) x - Inf)

    # bounds - if change things will need updated - manually
    lwr$init_log_Rt = rep(-10, length(lwr$init_log_Rt))
    lwr$log_Rt = rep(-10, n)
    lwr$log_Ft = rep(-15, n)
    lwr$F35 = -4.605
    lwr$F40 = -4.605
    lwr$F50 = -4.605
    lwr$sigmaR = 0.3
    lwr = lwr[names(p)] # make sure order is the same as pars
    lwr = unlist(lwr)

    upr$init_log_Rt = rep(10, length(upr$init_log_Rt))
    upr$log_Rt = rep(10, n)
    upr$log_Ft = rep(15, n)
    upr$F35 = 0
    upr$F40 = 0
    upr$F50 = 0
    upr$sigmaR = 10
    upr = upr[names(p)] # make sure order is the same as pars
    upr = unlist(upr)
    pars = p

    obj = RTMB::MakeADFun(f1, 
                         pars,
                         map = map)  

    fit = nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 control = list(iter.max=100000,
                                  eval.max=20000),
                 lower = lwr,
                 upper = upr)
    sdrep <- sdreport(obj, getJointPrecision = TRUE)

    list(obj = obj, sd = sdrep)

}

for(i in 1:10) {
    out = peels(model=f1, data=data, pars=pars2, map = list(sigmaR = factor(NA)), peel = i)
    saveRDS(out, file = here::here('retro', 'results', paste0('retro_' , i, '.rds')))
} 
peels(model=f1, data=data, pars=pars2, map = list(sigmaR = factor(NA)), peel = 1)
here::here()
peel_pars <- function(pars, peel) {
    p = pars2
    n = length(p$log_Rt)-peel
    p$log_Rt = head(p$log_Rt, -peel)
    p$log_Ft = head(p$log_Ft, -peel)
    upr = lapply(p, function(x) x + Inf)
    lwr = lapply(p, function(x) x - Inf)


    lwr$init_log_Rt = rep(-10, length(lwr$init_log_Rt))
    lwr$log_Rt = rep(-10, n)
    lwr$log_Ft = rep(-15, n)
    lwr$F35 = -4.605
    lwr$F40 = -4.605
    lwr$F50 = -4.605
    lwr$sigmaR = 0.3
    lwr = lwr[names(p)] # make sure order is the same as pars
    lwr = unlist(lwr)

    upr$init_log_Rt = rep(10, length(upr$init_log_Rt))
    upr$log_Rt = rep(10, n)
    upr$log_Ft = rep(15, n)
    upr$F35 = 0
    upr$F40 = 0
    upr$F50 = 0
    upr$sigmaR = 10
    upr = upr[names(p)] # make sure order is the same as pars
    upr = unlist(upr)
    list(pars = p, upr = upr, lwr = lwr)
    
}

peel_map <- function(map, pars){
  ## tricky part is map elements may or may not be specified so
  ## do this generically
  m <- map
  p <- pars
  stopifnot(is.list(map))
  for(i in names(m)){
    #print(i)
    #if(i=='log_q2_dev') browser()
    m[[i]] <- m[[i]][1:length(p[[i]])]
  }
  return(m)
}

run_retro <- function(year, folder, model, data, pars, peels = 10) {
    dir.create(here::here(year, folder, 'retro', 'results'), recursive=TRUE)

    # change model to run on retro data
    mdl = deparse(model)
    mdl = gsub("pars, data", "pars, r_data", mdl) 
    mdl = eval(parse(text = paste(mdl, collapse = "\n")))

    # define year to change
    for(i in 1:peels){
    # rename data 
    data = peel_data(data, peel)
    p = peel_pars(pars2, peel)
    pars = p$pars
    r_lwr = p$lwr
    r_upr = p$upr
    r_map = peel_map(map, pars)

    obj = RTMB::MakeADFun(f1, 
                         pars,
                         map = list(sigmaR = factor(NA)))  

    fit = nlminb(start = obj$par,
                   objective = obj$fn,
                   gradient = obj$gr,
                   control = list(iter.max=100000,
                                  eval.max=20000),
                   lower = lwr,
                   upper = upr)
    sdrep <- sdreport(obj, getJointPrecision = TRUE)

    

}
  return(list(obj = obj, fit = fit, sd = sdrep))
}

run_retro(year=2024, folder='retro', model=f1, data = data, pars = pars2, peels = 1)


    ifelse(data$fish_age_ind[]>(peel_n-1), 0, data$fish_age_ind) # never have comps the same year
    data[grep('catch_obs', names(data), value = TRUE)]
    

r_data$fish_age_ind[(peel_n-1):peel_n] <- 0
length(fish_age_ind)
data$fish_age_ind
    data %>% 
     filter(map(., ~ .xstarts_with('fish_age')))
    # change model to run on retro data
    mdl = deparse(f1)
    mdl = gsub("pars, data", "pars, r_data", mdl) 
    mdl = eval(parse(text = paste(mdl, collapse = "\n")))

    # drop some data
    i = 1
   r_data$years = r_data$years[r_data$years <= (max(r_data$years) - i)]

r_data$years <= (max(years - i))
}

fit_peel <- function(fit, peel, getsd=TRUE, ...) {
    obj <- RTMB::MakeADFun(mdl, 
                        pars,
                        map = map)  
}


retro_data = data
ff = deparse(f1)
retro_f = gsub("pars, data", "pars, retro_data", ff) 
head(retro_f)
deparse(retro_f)


grep('getAll(pars, data)', f1)
obj2$obj$env$data
obj2$env$data(post)
fit2$par
obj2$env
obj2$report(post)
obj$report(post[1,-ncol(post)])$Nat 
str(post)
fit2$data
obj2$env$data
ls(obj2$env)
obj2$par
obj$env$env
fit2$obj2$env$data
ls(fit2)
fit2$par
remotes::install_github("https://github.com/kaskr/RTMB", subdir="RTMB", force=TRUE)
1
