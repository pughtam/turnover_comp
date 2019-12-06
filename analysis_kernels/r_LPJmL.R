#========================================================================================
# This script calcualtes the forest turnover times for LPJmL
#----------------------------------------------------------------------------------------

# Require netcdf library
#----------------------------------------------------------------------------------------
if (!exists ('nc_open',   mode = 'function')) library (ncdf4)

# Read forest mask
#----------------------------------------------------------------------------------------
LPJmL_forest_mask    <- array (NA, dim = c (720, 360))
nc_name <- 'lpjml_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
LPJmL_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
LPJmL_forest_mask <- cbind (matrix (data = NA, ncol = 67, nrow = 720), LPJmL_forest_mask, matrix (data = NA, ncol = 14, nrow = 720)) # Add lines to get full North South extent
LPJmL_forest_mask [LPJmL_forest_mask == 2] <- NA # Set non-forest to NA for ease of analysis

# Read phenology mask
#----------------------------------------------------------------------------------------
LPJmL_pheno_mask  <- array (NA, dim = c (720, 360))
nc_name <- 'lpjml_cruncep_lai_annual_1901_2014_phenology_mask.nc'
ncin <-  nc_open (nc_name)
LPJmL_pheno_mask <- ncvar_get (ncin, 'phen_max_lai_phen_number')
LPJmL_pheno_mask <- cbind (matrix (data = NA, ncol = 67, nrow = 720), LPJmL_pheno_mask, matrix (data = NA, ncol = 14, nrow = 720)) # Add lines to get full North South extent

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------
LPJmL_CRUN_Cveg <- array (NA, dim = c (720, 360))
nc_name  <- 'lpjml_cruncep_cveg_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'cveg', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract only the years 1982 to 2011
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) 
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T)  
tmp.array3 <- cbind (matrix (data = NA, ncol = 67, nrow = 720), tmp.array2, matrix (data = NA, ncol = 14, nrow = 720)) # Add lines to get full North South extent
tmp.array3 [is.na (LPJmL_forest_mask)] <- NA
LPJmL_CRUN_Cveg <- tmp.array3

# Read CRUN NPP layer
#----------------------------------------------------------------------------------------
LPJmL_CRUN_NPP <- array (NA, dim = c (720, 360))
nc_name  <- 'lpjml_cruncep_npp_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'npp', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract only the years 1982 to 2011
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the 
tmp.array3 <- cbind (matrix (data = NA, ncol = 67, nrow = 720), tmp.array2, matrix (data = NA, ncol = 14, nrow = 720)) # Add lines to get full North South extent
tmp.array3 [is.na (LPJmL_forest_mask)] <- NA
LPJmL_CRUN_NPP <- tmp.array3 * 60.0 * 60.0 * 24.0 * 365.233333 # 7 /30 years in that period are leap years 

# Calculate CRUN tau_LPJmL
#----------------------------------------------------------------------------------------
LPJmL_CRUN_nta <- LPJmL_CRUN_Cveg / LPJmL_CRUN_NPP
LPJmL_CRUN_nta [is.na (LPJmL_forest_mask) | LPJmL_forest_mask == 2] <- NA

# Clean up
#----------------------------------------------------------------------------------------
rm (tmp.array, tmp.array1, tmp.array2)
#========================================================================================