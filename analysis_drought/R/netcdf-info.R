if (!require(RNetCDF))
  warning("package 'RNetCDF' is required to make use of these functions!")

#' return the dimension names and values in a list
#'
#' return the dimension names and values in a list
#' @param ncf either a open RNetCDF object from open.nc or a NetCDF file name.
#' @return list with 'names' and named values.
#' @author Joerg Steinkamp \email{joerg.steinkamp@@senckenberg.de}
#' @import RNetCDF
#' @export
dimensions.get.nc <- function(ncf=NULL) {
  if (class(ncf)=="character") {
    ncf <- open.nc(ncf)
  } else if (class(ncf)!="NetCDF") {
    stop(paste("Can't handle class'", class(ncf), "'!", sep=""))
  }
  #ret <- NULL

  dim.names <- c()
  dim.len <- list()
  i <- 0
  valid <- TRUE
  while(valid) {
    rval = tryCatch({
      dim.inq.nc(ncf, i)
    }, error = function(e) {
      return(FALSE)
    })
    if (is.logical(rval)) {
      valid=FALSE
    } else {
      dim.names <- append(dim.names, rval$name)
      eval(parse(text=paste("dim.len[['",rval$name,"']]=rval$length", sep="")))
    }
    i <- i+1
  }

  vars <- varnames.get.nc(ncf)

  ret <- list(names=dim.names)
  for (d in dim.names) {
    ## check if the dimension also exists as variable to get the values
    if (d %in% vars$name) {
      ret[[d]] = var.get.nc(ncf, d)
    } else {
      ret[[d]] = 1:dim.len[[d]]
    }
  }
  return(ret)
}

#' return the variable names and dimension ordering.
#'
#' return the variable names and dimension ordering.
#' @param ncf either a open RNetCDF object from open.nc or a NetCDF file name.
#' @param dim.names either NULL (default), FALSE or a vector of dimension names
#' @return list with variable 'names' in the NetCDF file. Either including dimension variables of not.
#' @author Joerg Steinkamp \email{joerg.steinkamp@@senckenberg.de}
#' @import RNetCDF
#' @export
varnames.get.nc <- function(ncf=NULL, dim.names=NULL) {
  if (class(ncf)=="character") {
    ncf <- open.nc(ncf)
  } else if (class(ncf)!="NetCDF") {
    stop(paste("Can't handle class'", class(ncf), "'!", sep=""))
  }

  if (is.logical(dim.names))
    if (!dim.names)
      dim.names <- dimensions.get.nc(ncf)$names

  var.names <- c()
  var.dims <- list()
  i <- 0
  valid <- TRUE
  while(valid) {
    rval = tryCatch({
      var.inq.nc(ncf, i)
    }, error = function(e) {
      return(FALSE)
    })
    if (is.logical(rval)) {
      valid=FALSE
    } else {
      var.names <- append(var.names, rval$name)
      eval(parse(text=paste("var.dims[['",rval$name,"']]=rval$dimids", sep="")))
    }
    i <- i+1
  }
  var.names <- setdiff(var.names, dim.names)

  for (d in dim.names) {
    var.dims[[d]] <- NULL
  }
  ret <- list(names=var.names, dims=var.dims)
  return(ret)
}

#' return NetCDF attributes.
#'
#' Either return global or variable specific attributes.
#'
#' @param ncf either a open RNetCDF object from open.nc or a NetCDF file name.
#' @param var.names either NULL (default) or a variable name.
#' @return list of variable attributes with id, name, type, length and value.
#' @author Joerg Steinkamp \email{joerg.steinkamp@@senckenberg.de}
#' @import RNetCDF
#' @export
attributes.get.nc <- function(ncf=NULL, var.name=NULL) {
  if (class(ncf)=="character") {
    ncf <- open.nc(ncf)
  } else if (class(ncf)!="NetCDF") {
    stop(paste("Can't handle class'", class(ncf), "'!", sep=""))
  }

  if (is.null(var.name))
    var.name <- "NC_GLOBAL"
  
  ret <- list()
  i <- 0
  valid <- TRUE
  while(valid) {
    rval = tryCatch({
      att.inq.nc(ncf, var.name, i)
    }, error = function(e) {
      return(FALSE)
    })
    if (is.logical(rval)) {
      valid=FALSE
    } else {
      ret[[i+1]] <- rval
      ret[[i+1]]$value = att.get.nc(ncf, var.name, ret[[i+1]]$name)
    }
    i <- i+1
  }
  return(ret)
}
