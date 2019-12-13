#========================================================================================
# This script calcualtes the forest turnover times for SEIB
#----------------------------------------------------------------------------------------

# Require netcdf library
#----------------------------------------------------------------------------------------
if (!exists ('nc_open',   mode = 'function')) library (ncdf4)

# Load maps library
#----------------------------------------------------------------------------------------
if (!exists ('forest_mask')) source ('R_scripts/read_and_plot_observational_data.R')

# Read forest mask
#----------------------------------------------------------------------------------------
SEIB_forest_mask    <- array (NA, dim = c (720, 360))
nc_name <- 'seib_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
SEIB_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
SEIB_forest_mask [SEIB_forest_mask == 2] <- NA
SEIB_forest_mask <- SEIB_forest_mask [,length (SEIB_forest_mask [1, ]):1]

# Read phenology mask
#----------------------------------------------------------------------------------------
SEIB_pheno_mask  <- array (NA, dim = c (720, 360))
nc_name <- 'seib_cruncep_lai_annual_1901_2014_phenology_mask.nc'
ncin <-  nc_open (nc_name)
SEIB_pheno_mask <- ncvar_get (ncin, 'phen_max_lai_phen_number')
SEIB_pheno_mask <- SEIB_pheno_mask [, length (SEIB_pheno_mask [1, ]):1]

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------
SEIB_CRUN_Cveg <- array (NA, dim = c (720, 360))
nc_name  <- 'seib_cruncep_cveg_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'cveg', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract only the years 1985 to 2014
tmp.array1 <- apply (tmp.array, c (1, 2, 3), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over time 
tmp.array3 <- tmp.array2 [, length (tmp.array2 [1, ]):1] # Turn the world up-side-down
tmp.array3 [is.na (forest_mask_NA) | is.na (SEIB_forest_mask)] <- NA # Apply the forest mask
SEIB_CRUN_Cveg <- tmp.array3

# Read CRUN NPP layer
#----------------------------------------------------------------------------------------
SEIB_CRUN_NPP <- array (NA, dim = c (720, 360))
nc_name  <- 'seib_cruncep_npp_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'npp', start = c (1, 1, 85), count = c (720, 360, 30)) # Extract only the years 1985 to 2014
tmp.array1 <- apply (tmp.array, c (1, 2), mean, na.rm = T) # Calculate mean over the 
tmp.array2 <- tmp.array1 [, length (tmp.array1 [1, ]):1] # Turn the world up-side-down
tmp.array2 [is.na (forest_mask_NA) | is.na (SEIB_forest_mask)] <- NA
SEIB_CRUN_NPP <- tmp.array2 * 60 * 60 * 24 * 365.23333

# Calculate CRUN tau SEIB
#----------------------------------------------------------------------------------------
SEIB_CRUN_nta <- SEIB_CRUN_Cveg / SEIB_CRUN_NPP
SEIB_CRUN_nta [is.na (SEIB_forest_mask) | SEIB_forest_mask == 2] <- NA

# Clean up
#----------------------------------------------------------------------------------------
rm (tmp.array, tmp.array1, tmp.array2)
#========================================================================================