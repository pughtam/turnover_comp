library(methods)
library(RNetCDF)

source("R/netcdf-info.R")

load("data/dm_events.RData")

time.lag <- 1

## LPj-wsl needs some reformating/conversation of the input data

## class for holding the observed and modeled data of a drought/mortality event (dme)
## as well as the wilcox significance test result
setClass("MortEvent",
         slots=c(id="character",
                 reference="character",
                 period="matrix",
                 gridlist="matrix",
                 ts="list",
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

## directory holding the mode outputs in with the subdirectories <model>/<NetCDF files>
nc.dir <- "/senckenberg.de/cub/more_akira/Tom_Pugh_Drought_Mortality/data"

## hierachically organized list of model name with sub-lists of climate drivers and
## and file pattern. Those files are summed up.
models <- list("CABLE-POP"=list(drivers=c("CRUNCEP", "IPSL-CM5A-LR"),
                                file.patterns=c("cmortwooddist",
                                                "cmortwoodresourcecrowding",
                                                "cmortwoodresourcelim")),
               "JULESC2"=list(drivers=c("CRUNCEP", "IPSL"),
                              file.patterns=c("wood_litc",
                                              "lit_c_dyn_wood",
                                              "lit_c_dyn_leaf",
                                              "lit_c_dyn_root")),
               "LPJ-GUESS"=list(drivers=c("CRUNCEP", "IPSL"),
                                file.patterns=c("cmort_badallom",
                                                "cmort_bioclim",
                                                "cmort_dist",
                                                "cmort_fire",
                                                "cmort_groweff",
                                                "cmort_mmin",
                                                "cmort_nbiom")),
               "LPJmL"=list(drivers=c("CRUNCEP", "IPSL-CM5A_LR"),
                            file.patterns=c("cmort_fire",
                                            "cmort_growth_efficiency",
                                            "cmort_heatstress",
                                            "cmort_neg_allocation",
                                            "cmort_neg_biomass",
                                            "cmort_shade")),
               "LPJ-wsl"=list(drivers=c("CRUNCEP", "IPSL-CM5A-LR"),
                              file.patterns=c("cmort_fire_tree",
                                              "cmort_greff_tree",
                                              "cmort_bioclim_tree",
                                              "cmort_heat_tree",
                                              "cmort_space_tree",
                                              "cmort_fire_grass",
                                              "cmort_greff_grass",
                                              "cmort_bioclim_grass",
                                              "cmort_heat_grass",
                                              "cmort_space_grass")),
               "ORCHIDEE"=list(drivers=c("CRUNCEP", "IPSL-CM5A-LR"),
                               file.patterns=c("total_turn")),
               "SEIB-DGVM"=list(drivers=c("CRUNCEP", "ipsl"),
                                file.patterns=c("cmort_etc",
                                                "cmort_fire_atmosp",
                                                "cmort_fire_litter",
                                                "cmort_gap",
                                                "cmort_greff",
                                                "cmort_heat",
                                                "cmort_lim")))
model.summary <- list()
model.results <- list()
## set a random seed for reproducability
set.seed(987654)
for (m in names(models)) {
  message(paste("***", m, "***"))
  in.dir <- file.path(nc.dir, m)
  for (d in models[[m]][["drivers"]]) {
    message(paste0(d, " :"), appendLF = FALSE)

    ## get the forest mask
    if (m == "CABLE-POP") {
      fname <- file.path(nc.dir, paste0(tolower(m), "_cruncep_lai_annual_1901_2015_forest_mask_v4.nc"))
    } else if (m == "JULESC2") {
      fname <- file.path(nc.dir, paste0(tolower(m), "_cruncep_lai_annual_1901_2014_forest_mask_v4_origin_05.nc"))
    } else {
      fname <- file.path(nc.dir, paste0(tolower(m), "_cruncep_lai_annual_1901_2014_forest_mask_v4.nc"))
    }
    ncin.mask <- open.nc(fname)
    dims.mask <- dimensions.get.nc(ncin.mask)
    
    files <- list.files(in.dir, pattern=paste0(tolower(d), "_cveg_", c("annual", "year") , collapse="|"), full.names=TRUE)
    ## some additional filtering
    files = files[!grepl("gz$", files)]
    fname = files[!grepl("seib_cruncep_cveg_annual_1901_2099.nc4", files)]

    if (length(fname) != 1)
      warning("*** more than one cveg file!")

    ncin.cveg <- open.nc(fname)
    dims.cveg <- dimensions.get.nc(ncin.cveg)
    lon.name  <- dims.cveg$name[tolower(dims.cveg$name) %in% c("lon", "longitude", "x")]
    lat.name  <- dims.cveg$name[tolower(dims.cveg$name) %in% c("lat", "latitude", "y")]
    time.name <- dims.cveg$name[tolower(dims.cveg$name) %in% c("time_counter", "time", "t", "z")]

    ## check agreement of data and mask grid-layout
    if (!all(dims.cveg[[lon.name]] == dims.mask[[lon.name]])) {
      stop("Longitude of cveg and mask do not match!")
    }
    if (!all(dims.cveg[[lat.name]] == dims.mask[[lat.name]])) {
      if (all(dims.cveg[[lat.name]] == rev(dims.mask[[lat.name]]))) {
        stop("Latitude of cveg and mask inverted!")
      } else {
        stop("Latitude of cveg and mask do not match!")
      }
    }
    
    if (tolower(d) == "cruncep" && m != "CABLE-POP") {
      message(paste0("shift one gridcell north (",m, "/",d,")"))
      dm.grid.ids <- lapply(dm.events, find.grid.ids, lon=dims.cveg[[lon.name]], lat=dims.cveg[[lat.name]] - 0.5)
    } else {
      dm.grid.ids <- lapply(dm.events, find.grid.ids, lon=dims.cveg[[lon.name]], lat=dims.cveg[[lat.name]])
    }
    
    mortTS <- list()
    ## loop over individual files
    for (p in models[[m]]$file.patterns) {
      fname <- list.files(in.dir, pattern=paste0(tolower(d), "_", p, "_"), full.names=TRUE)
      fname = fname[!grepl("gz$", fname)]
      if (length(fname) != 1)
        warning("*** more than one mort file!")

      message(paste0(" ", p, " "), appendLF=FALSE)
      ncin.mort <- open.nc(fname)
      vars <- varnames.get.nc(ncin.mort)
      time <- var.get.nc(ncin.mort, time.name)
      monthly <- FALSE
      if (length(time) > 1000)
        monthly <- TRUE

      mortTS[[p]] <- lapply(names(dm.events), function(x) {
        ## mortTS[[p]] = foreach (x=names(dm.events)) %dopar% {
        message(x)
        ## get the rectangular data
        mask <- var.get.nc(ncin.mask, "forest_30yr_any_10_years",
                           start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1]),
                           count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2]))
        mask[!is.finite(mask)] = 3
        if (length(vars$dims[[p]]) == 3) {
          cveg.vals <- var.get.nc(ncin.cveg, "cveg",
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA))
          mort.vals <- var.get.nc(ncin.mort, p,
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA))
          ## eliminate missing values in LPJ-wsl
          mort.vals[abs(mort.vals + 99999/365/86400) < 1.e-18] = NA
          ## set negative mortality fluxes to zero
          mort.vals[mort.vals < 0] = 0.0
        } else if (length(vars$dims[[p]]) == 4) {
          cveg.vals <- var.get.nc(ncin.cveg, "cveg",
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA, NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA, NA))
          mort.vals <- var.get.nc(ncin.mort, p,
                                  start=c(dm.grid.ids[[x]][1,1], dm.grid.ids[[x]][2,1], NA, NA),
                                  count=c(dm.grid.ids[[x]][1,2], dm.grid.ids[[x]][2,2], NA, NA))
          ## PFT is the second last dimension for all models without SEIB-DGVM, here it is the last
          ## (summing over it)
          if (m=="SEIB-DGVM") {
            pftDim <- length(dim(mort.vals))
          } else {
            pftDim <- length(dim(mort.vals)) - 1
          }
          otherDim <- 1:length(dim(mort.vals))
          otherDim <- otherDim[otherDim != pftDim]
          ## set negative mortality fluxes to zero
          mort.vals[mort.vals < 0] = 0.0
          cveg.vals = apply(cveg.vals, otherDim, sum)
          mort.vals = apply(mort.vals, otherDim, sum)
        } else {
          stop("Don't know how to handle data with dimensions other than 3 or 4!")
        }
        ## apply the forest mask and the mask of the gridpoints within the rectangle
        lon <- dims.cveg[[lon.name]][dm.grid.ids[[x]][1,1]:(dm.grid.ids[[x]][1,1] + dm.grid.ids[[x]][1,2] - 1)]
        lat <- dims.cveg[[lat.name]][dm.grid.ids[[x]][2,1]:(dm.grid.ids[[x]][2,1] + dm.grid.ids[[x]][2,2] - 1)]
        mask <- array(mask, c(length(lon), length(lat)))
        rect.grid <- data.frame(lon=rep(lon, length(lat)), lat=rep(lat, each=length(lon)), mask=as.vector(mask))
        for (i in 1:nrow(rect.grid)) {
          for (j in 1:ncol(dm.events[[x]]$gridlist)) {
            if (rect.grid$mask[i] == 1 &&
                rect.grid$lon[i] == dm.events[[x]]$gridlist[1,j] &&
                rect.grid$lat[i] == dm.events[[x]]$gridlist[2,j]) {
              rect.grid$mask[i] = TRUE
            } else {
              rect.grid$mask[i] = FALSE
            }
          }
        }

        ## simple mean, not area weighted. Is ok, as long as the lat extent isn't too big.
        if (length(dim(mort.vals)) > 1) {
          cveg.vals = apply(cveg.vals, length(dim(cveg.vals)), function(y){mean(y[rect.grid$mask], na.rm=TRUE)})
          mort.vals = apply(mort.vals, length(dim(mort.vals)), function(y){mean(y[rect.grid$mask], na.rm=TRUE)})
        } else if (mask != 1) {
          mort.vals[] = NA
        }
        if (monthly) {
          year.id = rep(1:(length(time) / 12), each=12)
          mort.vals = aggregate(mort.vals, list(year.id), weighted.mean, w=c(31,28,31,30,31,30,31,31,30,31,30,31), na.rm=TRUE)[, 2]
        }
        
        ## introduce a cveg limit threshold
        values <- cveg.vals
        ## values[cveg.vals < 1.0] = NA
        values <- mort.vals * 365.0 * 86400.0 / values
        values[!is.finite(values)] = NA
        ### HARDCODED number of years and start year!!!
        return(ts(data.frame(mort.frac=values[1:114], mort.flux=mort.vals[1:114], cveg=cveg.vals[1:114]), start=1901, freq=1))
      }) ## end lapply over dm.events
      names(mortTS[[p]]) <- names(dm.events)
    } ## loop over mortality files

    ## message("mortality timeseries extracted")
    
    mortEvent <- lapply(names(dm.events), function(x) {
      tmp <- new("MortEvent", id=x,
                 reference=dm.events[[x]][["reference"]],
                 period=dm.events[[x]][["period"]],
                 gridlist=dm.events[[x]][["gridlist"]])
      tmp@ts=list()
      for (i in 1:length(models[[m]]$file.patterns)) {
        if (i == 1) {
          tmp@ts[["Total"]] = mortTS[[models[[m]]$file.patterns[i]]][[x]]
        } else {
          tmp@ts[["Total"]][,1:2] = tmp@ts[["Total"]][,1:2] + mortTS[[models[[m]]$file.patterns[i]]][[x]][, 1:2]
        }
        tmp@ts[[models[[m]]$file.patterns[i]]] = mortTS[[models[[m]]$file.patterns[i]]][[x]]
      }
      colnames(tmp@ts[["Total"]]) <- colnames(tmp@ts[[models[[m]]$file.patterns[1]]])
      return(tmp)
    })
    names(mortEvent) <- names(dm.events)

    ## message("mortEvents created")
    
    ## problem with zeroes. According to Bob add a tiny amount to the zeroes and log-transform
    ## adding 5 subsequent years
    ## the wilcox test is insensitive to the position of the zero transformation
    for (dme in names(mortEvent)) {
      ##message(dme)
      mortEvent[[dme]]@wilcox = sapply(names(mortEvent[[dme]]@ts), function(x) {
        dme.ts <- mortEvent[[dme]]@ts[[x]][,'mort.frac']
        dme.ts[dme.ts < 0] = NA
        if (length(unique(dme.ts)) < 2)
          return(NA)
        if (sum(is.finite(dme.ts)) < 2)
          return(NA)
        if (min(dme.ts, na.rm=TRUE) == 0) {
          delta = min(dme.ts[dme.ts > 0], na.rm=TRUE) / 5
          full = dme.ts + delta
        } else {
          full = dme.ts
        }
        sample = full[time(full) >= mortEvent[[dme]]@period[1] + time.lag & time(full) <= mortEvent[[dme]]@period[2] + time.lag]
        if (sum(is.finite(sample)) < 2)
          return(NA)
        full = as.vector(full)
        wt = wilcox.test(log10(full), log10(sample), alternative="less")
        return(wt$p.value)
      })

      ## message("random")
      
      ## return as list
      mortEvent[[dme]]@random.wilcox = lapply(names(mortEvent[[dme]]@ts), function(x) {
        dme.ts <- mortEvent[[dme]]@ts[[x]][,'mort.frac']
        dme.ts[dme.ts < 0] = NA
        if (length(unique(dme.ts)) < 2)
          return(NA)
        if (sum(is.finite(dme.ts)) < 2)
          return(NA)
        res = NULL
        if (min(dme.ts, na.rm=TRUE) == 0) {
          delta = min(dme.ts[dme.ts > 0], na.rm=TRUE) / 5
          full = dme.ts + delta
        } else {
          full = dme.ts
        }
        years = min(time(full)):(max(time(full)) - (mortEvent[[dme]]@period[2] - mortEvent[[dme]]@period[1]))
        for (y in sample(years, 10, replace=TRUE)) {
          sample = full[time(full) > y & time(full) <= y + (mortEvent[[dme]]@period[2] - mortEvent[[dme]]@period[1])]
          if (sum(is.finite(sample)) < 2) {
            res = append(res, NA)
          } else {
            wt = wilcox.test(log10(as.vector(full)), log10(sample), alternative="less")
            res = append(res, wt$p.value)
          }
        }
        return(res)
      })
      names(mortEvent[[dme]]@random.wilcox) <- names(mortEvent[[dme]]@ts)
    }

    ## extract the wilcox test for true and random time period
    ## per mortality process
    for (p in names(mortEvent[[1]]@ts)) {
      match = NULL
      random.match = NULL
      for (dme in names(dm.events)) {
        if (length(unique(mortEvent[[dme]]@ts[[p]][,'mort.frac'])) < 2) {
          next
        }
        match = append(match, mortEvent[[dme]]@wilcox[p])
        random.match = append(random.match, mortEvent[[dme]]@random.wilcox[[p]])
      }

      model.summary[[m]][[d]] <- rbind(model.summary[[m]][[d]], c(n.match=sum(match < 0.05, na.rm=TRUE),
                                                                  n.total=length(match),
                                                                  event=sum(match < 0.05, na.rm=TRUE) / length(match),
                                                                  n.match.random=sum(random.match < 0.05, na.rm=TRUE),
                                                                  n.total.random=length(random.match),
                                                                  random=sum(random.match < 0.05, na.rm=TRUE) / length(random.match)))
    }
    row.names(model.summary[[m]][[d]]) <- names(mortEvent[[1]]@ts)
    model.results[[m]][[d]] = mortEvent
    message("done.")
  } ## end loop over drivers
} ## end loop over DGVMs


## renaming driver IPSL to be equal for all models
for (m in names(model.results)) {
  dnames <- names(model.results[[m]])
  dnames[grepl("ipsl", tolower(dnames))] = "IPSL-CM5A-LR"
  names(model.results[[m]]) <- dnames
}
rm(dnames)

## remove the prefix "cmort" 
for (m in names(model.results)) {
  for (d in names(model.results[[m]])) {
    for (l in names(model.results[[m]][[d]])) {
      if (m == 'CABLE-POP') {
        pnames <- names(model.results[[m]][[d]][[l]]@ts)
        pnames[pnames == "cmort"] = "cmort_cmort"
        pnames = sub("cmort|cmort_", "", pnames)
        names(model.results[[m]][[d]][[l]]@ts) <- pnames

        pnames <- attr(model.results[[m]][[d]][[l]]@wilcox, 'name')
        pnames[pnames=="cmort"] = "cmort_cmort"
        pnames = sub("cmort|cmort_", "", pnames)
        attr(model.results[[m]][[d]][[l]]@wilcox, 'names') <- pnames

        pnames <- names(model.results[[m]][[d]][[l]]@random.wilcox)
        pnames[pnames=="cmort"] = "cmort_cmort"
        pnames = sub("cmort|cmort_", "", pnames)
        names(model.results[[m]][[d]][[l]]@random.wilcox) <- pnames
      } else {
        pnames <- names(model.results[[m]][[d]][[l]]@ts)
        pnames = sub("cmort|cmort_", "", pnames)
        names(model.results[[m]][[d]][[l]]@ts) <- pnames

        pnames <- attr(model.results[[m]][[d]][[l]]@wilcox, 'name')
        pnames = sub("cmort|cmort_", "", pnames)
        attr(model.results[[m]][[d]][[l]]@wilcox, "names") <- pnames

        pnames <- names(model.results[[m]][[d]][[l]]@random.wilcox)
        pnames = sub("cmort|cmort_", "", pnames)
        names(model.results[[m]][[d]][[l]]@random.wilcox) <- pnames
      }
    }  
  }
}
rm(pnames)

save(file="data/mortEvents.RData", compress="xz", model.summary, model.results)
