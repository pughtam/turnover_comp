#========================================================================================
# This script calcualtes the forest turnover times for JULES.
#----------------------------------------------------------------------------------------

# Read forest mask
#----------------------------------------------------------------------------------------
JULES_forest_mask    <- array (NA, dim = c (192, 112+29+4))
nc_name <- 'Shared_files/forest_mask/JULESC2_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
tmp.array1 <- ncvar_get (ncin, 'forest_30yr_any_10_years')
tmp.array2 <- cbind (matrix (rep (NA, 192 * 29), nrow = 192), rbind (tmp.array1 [97:192, ], tmp.array1 [1:96, ]), matrix (rep (NA, 192 * 4), nrow = 192))
JULES_forest_mask <- tmp.array2
JULES_forest_mask [JULES_forest_mask == 2] <- NA

# Read phenology mask
#----------------------------------------------------------------------------------------
JULES_pheno_mask  <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/phenology/JULESC2_cruncep_lai_annual_1901_2014_phenology_mask.nc'
ncin <-  nc_open (nc_name)
tmp.array1       <- ncvar_get (ncin, 'phen_max_lai_phen_number')
JULES_pheno_mask <- cbind (matrix (rep (NA, 192 * 29), nrow = 192), rbind (tmp.array1 [97:192, ], tmp.array1 [1:96, ]), matrix (rep (NA, 192 * 4), nrow = 192))

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------
JULES_CRUN_Cveg <- array (NA, dim = c (192, 112+29+4))
nc_name  <- 'data/JULES/cruncep/JULESC2_CRUNCEP_cveg_gb_Annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'cveg_gb', start =  c (1, 1, 85), count = c (192, 112, 30)) # Extract variable
tmp.array1 <- apply (tmp.array, c (1, 2), mean, na.rm = T) # Calculate mean over the 
tmp.array2 <- cbind (matrix (rep (NA, 192 * 29), nrow = 192), rbind (tmp.array1 [97:192, ], tmp.array1 [1:96, ]), matrix (rep (NA, 192 * 4), nrow = 192))
JULES_CRUN_Cveg <- tmp.array2

# Read CRUN NPP layer
#----------------------------------------------------------------------------------------
JULES_CRUN_NPP <- array (NA, dim = c (192, 112+29+4))
nc_name  <- 'data/JULES/cruncep/JULESC2_CRUNCEP_npp_gb_Annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'npp_gb', start = c (1, 1, 85), count = c (192, 112, 30)) # Extract variable
tmp.array1 <- apply (tmp.array, c (1, 2), mean, na.rm = T) # Calculate mean over the 
tmp.array2 <- cbind (matrix (rep (NA, 192 * 29), nrow = 192), rbind (tmp.array1 [97:192, ], tmp.array1 [1:96, ]), matrix (rep (NA, 192 * 4), nrow = 192))
JULES_CRUN_NPP <- tmp.array2 * 60.0 * 60.0 * 24.0 * (365 + 7/30) # 7 /30 years in that period are leap years 

# Calculate CRUN tau_NPP for JULES
#----------------------------------------------------------------------------------------
JULES_CRUN_nta <- JULES_CRUN_Cveg / JULES_CRUN_NPP
JULES_CRUN_nta [is.na (JULES_forest_mask) | JULES_forest_mask == 2] <- NA

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------
JULES_CRUN_nta_adj <- JULES_CRUN_nta * (forest_mask_julesgrid_NA / 100.0) 

# Clean up
#----------------------------------------------------------------------------------------
rm (tmp.array, tmp.array1, tmp.array2)
#========================================================================================