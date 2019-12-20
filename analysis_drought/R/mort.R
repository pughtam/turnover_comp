library(methods)
library(ggplot2)
library(RNetCDF)
source("R/netcdf-info.R")

## class for holding the observed and modeled data of a drought/mortality event (dme)
## as well as the wilcox significance test result
setClass("MortEvent",
         slots=c(id="character",
                 reference="character",
                 period="matrix",
                 gridlist="matrix",
                 mort="list",
                 wilcox="vector",
                 random.wilcox="vector"))

## find the position of the rectangular extent in the drought/mortality events list
find.grid.ids <- function(x, lon, lat) {
  xmin <- min(which(lon==min(x$gridlist[1,])), which(lon==max(x$gridlist[1,])))
  xmax <- max(which(lon==min(x$gridlist[1,])), which(lon==max(x$gridlist[1,])))
  ymin <- min(which(lat==min(x$gridlist[2,])), which(lat==max(x$gridlist[2,])))
  ymax <- max(which(lat==min(x$gridlist[2,])), which(lat==max(x$gridlist[2,])))
  return(matrix(c(xmin, xmax - xmin + 1,
                  ymin, ymax - ymin + 1), ncol=2, byrow=TRUE))
}

## directory holding the mode outputs in with the subdirectories <model>/<driver>/<NetCDF files>
nc.dir <- "/senckenberg.de/cub/bigdata/DroughtMortality_Pugh/data"

## JULES data has no valid time axis and is coarser than 0.5 degree

## hierachically organized list of model name with sub-lists of climate drivers and
## and file pattern. Those files are summed up.
models <- list("LPJ-GUESS"=list(drivers=c("CRUNCEP", "IPSL"),
                                file.patterns=c("cmort_badallom", "cmort_bioclim",
                                                "cmort_dist", "cmort_fire",
                                                "cmort_groweff", "cmort_mmin", "cmort_nbiom")),
               "CABLE-POP"=list(drivers=c("CRUNCEP", "IPSL-CM5A-LR"),
                                file.patterns=c("cmort", "cmortfineroot", "cmortleaf", "cmortwooddist",
                                                "cmortwood", "cmortwoodresourcecrowding", "cmortwoodresourcelim")),
               ##               "JULESC2"=list(drivers=c("CRUNCEP", "IPSL"),
               ##                            file.patterns=c("cmort_comp","cmort_dist")),
               "LPJmL"=list(drivers=c("CRUNCEP", "IPSL-CM5A_LR"),
                            file.patterns=c("cmort_fire", "cmort_growth_efficiency", "cmort_heatstress",
                                            "cmort_neg_allocation", "cmort_neg_biomass",
                                            "cmort_redroduction_cost", "cmort_shade",
                                            "cmort_turnover_leaf", "cmort_turnover_root")),
               ##               "LPJ-wsl"=list(drivers="CRUNCEP",
               ##                              file.patterns=c("cmort_bioclim", "mort_greff",
               ##                                              "mort_space","mort_fire","mort_heat")),
               ##               "ORCHIDEE"=list(drivers="CRUNCEP",
               ##                                file.patterns=c("MORTALITY")),
               "SEIB-DGVM"=list(drivers=c("CRUNCEP", "ipsl"),
                                file.patterns=c("cmort_etc", "cmort_fire_atmosp",
                                                "cmort_fire_litter", "cmort_gap",
                                                "cmort_greff", "cmort_heat",
                                                "cmort_lim")))
model.summary <- list()
model.results <- list()
## set a random seed for reproducability
set.seed(98765)
for (m in names(models)) {
  message(m)
  in.dir <- file.path(nc.dir, m)[1]
  for (d in models[[m]][["drivers"]]) {
    message(d)
    load("data/dm_events.RData")
    files <- list.files(in.dir, pattern=paste0(tolower(d), "_cveg_", c("annual", "year") , collapse="|"), full.names=TRUE)
    ## some additional filtering
    files = files[!grepl("gz$", files)]
    fname = files[!grepl("seib_cruncep_cveg_annual_1901_2014.nc4", files)]

    if (length(fname) != 1)
      warning("*** more than one cveg file!")

    ncin.cveg <- open.nc(fname)
    dims      <- dimensions.get.nc(ncin.cveg)
    lon.name  <- dims$name[tolower(dims$name) %in% c("lon", "longitude", "x")]
    lat.name  <- dims$name[tolower(dims$name) %in% c("lat", "latitude", "y")]
    time.name <- dims$name[tolower(dims$name) %in% c("time_counter", "time", "t", "z")]

    dm.grid.ids <- lapply(dm.events, find.grid.ids, lon=dims[[lon.name]], lat=dims[[lat.name]])
    
    ## getting informations from the first mortality file
                                        #files <- list.files(in.dir, pattern=paste0(tolower(d), "_", models[[m]][["file.patterns"]], collapse="|"), full.names=TRUE)
                                        #files = files[!grepl("gz$", files)]
                                        #print(files)

    mortTS <- list()
    ## loop over individual files
    for (p in models[[m]]$file.patterns) {
      fname <- list.files(in.dir, pattern=paste0(tolower(d), "_", p, "_"), full.names=TRUE)
      fname = fname[!grepl("gz$", fname)]
      if (length(fname) != 1)
        warning("*** more than one cmort file!")

      ncin.mort <- open.nc(fname)
      vars <- varnames.get.nc(ncin.mort)
      time <- var.get.nc(ncin.mort, time.name)
      monthly <- FALSE
      if (length(time) > 1000)
        monthly <- TRUE

      ## get the lat/lon grid IDs per event as list


      mortTS[[p]] <- lapply(names(dm.events), function(x) {
        ## get the rectangular data
        if (length(vars$dims[[p]]) == 3) {
          cveg.vals <- var.get.nc(ncin.cveg, "cveg",
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA))
          mort.vals <- var.get.nc(ncin.mort, p,
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA))
        } else if (length(vars$dims[[p]]) == 4) {
          cveg.vals <- var.get.nc(ncin.cveg, "cveg",
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA, NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA, NA))
          mort.vals <- var.get.nc(ncin.mort, p,
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA, NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA, NA))
          ## print(paste("JS_DEBUG", x, paste0(dim(values), collapse=", ")))
          ## PFT is the second last dimension for all models without SEIB-DGVM, here it is the last
          ## (summing over it)
          if (m=="SEIB-DGVM") {
            pftDim <- length(dim(mort.vals))
          } else {
            pftDim <- length(dim(mort.vals)) - 1
          }
          otherDim <- 1:length(dim(mort.vals))
          otherDim <- otherDim[otherDim != pftDim]
          cveg.vals = apply(cveg.vals, otherDim, sum)
          mort.vals = apply(mort.vals, otherDim, sum)
        } else {
          stop("Don't know how to handle data with dimensions other than 3 or 4!")
        }
        ## apply a mask of the gridpoints within the rectangle
        lon <- dims[[lon.name]][dm.grid.ids[[x]][1,1]:(dm.grid.ids[[x]][1,1] + dm.grid.ids[[x]][1,2] - 1)]
        lat <- dims[[lat.name]][dm.grid.ids[[x]][2,1]:(dm.grid.ids[[x]][2,1] + dm.grid.ids[[x]][2,2] - 1)]
        rect.grid <- data.frame(lon=rep(lon, length(lat)), lat=rep(lat, each=length(lon)), mask=FALSE)
        for (i in 1:nrow(rect.grid)) {
          for (j in 1:ncol(dm.events[[x]]$gridlist)) {
            if (rect.grid$lon[i] == dm.events[[x]]$gridlist[1,j] &&
                rect.grid$lat[i] == dm.events[[x]]$gridlist[2,j]){
              rect.grid$mask[i]=TRUE
              next
            }
          }
        }
        ## simple mean, not area weighted. Is ok, as long as the lat extent isn't too big.
        if (length(dim(mort.vals)) > 1) {
          cveg.vals = apply(cveg.vals, length(dim(cveg.vals)), function(y){mean(y[rect.grid$mask], na.rm=TRUE)})
          mort.vals = apply(mort.vals, length(dim(mort.vals)), function(y){mean(y[rect.grid$mask], na.rm=TRUE)})
        }
        if (monthly) {
          year.id = rep(1:(length(time) / 12), each=12)
          mort.vals = aggregate(mort.vals, list(year.id), weighted.mean, w=c(31,28,31,30,31,30,31,31,30,31,30,31), na.rm=TRUE)[2]
        }
        return(mort.vals * 365 * 86400 / cveg.vals)
        
      }) ## end lapply over dm.events
      names(mortTS[[p]]) <- names(dm.events)
    }

    mortEvent <- lapply(names(dm.events), function(x) {
      tmp <- new("MortEvent", id=x, reference=dm.events[[x]][["reference"]],
                 period=dm.events[[x]][["period"]],
                 gridlist=dm.events[[x]][["gridlist"]])
      for (i in 1:length(models[[m]]$file.patterns)) {
        tmp@mort[[models[[m]]$file.patterns[i]]] = ts(mortTS[[models[[m]]$file.patterns[i]]][[x]], start=1901, freq=1)
        if (i == 1) {
          tmp@mort[["Total"]] = tmp@mort[[models[[m]]$file.patterns[i]]]
        } else {
          tmp@mort[["Total"]] = tmp@mort[["Total"]] + tmp@mort[[models[[m]]$file.patterns[i]]]
        }
      }
      return(tmp)
    })
    names(mortEvent) <- names(dm.events)

    ## problem with zeroes. According to Bob add a tiny amount to the zeroes and log-transform
    ## adding 5 subsequent years
    ## the wilcox test is insensitive to the position of the zero transformation
    for (dme in names(mortEvent)) {
      mortEvent[[dme]]@wilcox = sapply(names(mortEvent[[dme]]@mort), function(x) {
        dme.ts <- mortEvent[[dme]]@mort[[x]]
        if (length(unique(dme.ts)) < 2)
          return(NA)
        delta = min(dme.ts[dme.ts > 0]) / 5
        full = dme.ts + delta
        sample = full[time(full) > mortEvent[[dme]]@period[1] & time(full) <= mortEvent[[dme]]@period[2] + 5]
        full = as.vector(full)
        wt = wilcox.test(log10(full), log10(sample), alternative="less")
        return(wt$p.value)
      })

      ## return as list
      mortEvent[[dme]]@random.wilcox = lapply(names(mortEvent[[dme]]@mort), function(x) {
        dme.ts <- mortEvent[[dme]]@mort[[x]]
        if (length(unique(dme.ts)) < 2)
          return(NA)
        res = NULL
        delta = min(dme.ts[dme.ts > 0]) / 5
        full = dme.ts + delta
        years = min(time(full)):(max(time(full))-5)
        for (y in sample(years, 10, replace=TRUE)) {
          sample = full[time(full) > y & time(full) <= y + 5 + (mortEvent[[dme]]@period[2] - mortEvent[[dme]]@period[1])]
          wt = wilcox.test(log10(as.vector(full)), log10(sample), alternative="less")
          
          res = append(res, wt$p.value)
        }
        return(res)
      })
      names(mortEvent[[dme]]@random.wilcox) <- names(mortEvent[[dme]]@mort)
    }

    ## extract the wilcox test for true and random time period
    ## per mortality process
    for (p in names(mortEvent[[1]]@mort)) {
      match = NULL
      random.match = NULL
      for (dme in names(dm.events)) {
        if (length(unique(mortEvent[[dme]]@mort[[p]])) < 2) {
          next
        }
        match = append(match, mortEvent[[dme]]@wilcox[p])
        random.match = append(random.match, mortEvent[[dme]]@random.wilcox[[p]])
      }
      
      model.summary[[m]][[d]] <- rbind(model.summary[[m]][[d]], c(event=sum(match<0.05)/length(match),
                                                        random=sum(random.match<0.05)/length(random.match)))
    }
    row.names(model.summary[[m]][[d]]) <- names(mortEvent[[1]]@mort)
    model.results[[m]][[d]] = mortEvent
    message("done.")
    
  } ## end loop over drivers
} ## end loop over DGVMs

save.image(file="data/mortEvents.RData", compress="xz")
