library(xlsx)

## turn warnings into errors
options(warn=2)

str2num.lonlat <- function(str, warn.only=FALSE) {
  rv <- NULL
  for (i in 1:length(str)) {
    if (grepl('W$', str[i])) {
      rv <- append(rv, -1 * as.numeric(substr(str[i], 1, nchar(str[i])-1)))
    } else if (grepl('E$', str[i])) {
      rv <- append(rv, as.numeric(substr(str[i], 1, nchar(str[i])-1)))
    } else if (grepl('S$', str[i])) {
      rv <- append(rv, -1 * as.numeric(substr(str[i], 1, nchar(str[i])-1)))
    } else if (grepl('N$', str[i])) {
      rv <- append(rv, as.numeric(substr(str[i], 1, nchar(str[i])-1)))
    } else if (warn.only) {
      warning(paste0("No N/S or E/W suffix given. Assume real number ('", str[i],"')."))
      rv <- append(rv, as.numeric(str[i]))
    } else {
      stop(paste0("No N/S or E/W suffix given. Assume real number ('", str[i],"')."))
    }
  }
  return(rv)
}

get.all.drought.mortality.events <- function(raw.events) {
  events <- apply(raw.events, 1, function(x) {
    pos <- strsplit(x['GC'], ";")
    coords <- strsplit(pos[[1]], "/")
    coords <- sapply(coords, function(y) {
      if (grepl('[WE]$', y[1])) {
        lon <- str2num.lonlat(y[1])
      } else if (grepl('[WE]$', y[2])) {
        lon <- str2num.lonlat(y[2])
      }
      if (grepl('[NS]$', y[1])) {
        lat <- str2num.lonlat(y[1])
      } else if (grepl('[NS]$', y[2])) {
        lat <- str2num.lonlat(y[2])
      }

      return(matrix(c(lon, lat), ncol=2))
    })
    return(list(period=matrix(as.numeric(c(x['start'], x['end']))),
                gridlist=coords,
                reference=x['reference']))
  })
  names(events) <- raw.events$TblID
  return(events)
}

raw.events <- read.xlsx("data/locations_new_modSD.xlsx", "Data")
references <- read.xlsx("data/locations_new_modSD.xlsx", "References")
colnames(raw.events) <- c("TblID", "start", "end", "reference", "GC")
raw.events$TblID <- as.character(raw.events$TblID)
raw.events$reference <- as.character(raw.events$reference)
raw.events$GC <- as.character(raw.events$GC)
  
## remove rows with missing data (header + one multiple occurence)
raw.events <- raw.events[complete.cases(raw.events), ]
dm.events <- get.all.drought.mortality.events(raw.events)
save(file="data/dm_events.RData", dm.events, references)

