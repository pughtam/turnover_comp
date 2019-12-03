#========================================================================================
# This script calcualtes the forest turnover times for CABLE-POP.
#----------------------------------------------------------------------------------------

# Read forest mask
#----------------------------------------------------------------------------------------
CABLE_forest_mask <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/forest_mask/CABLE-POP_cruncep_lai_annual_1901_2015_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
CABLE_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
CABLE_forest_mask [CABLE_forest_mask == 2] <- NA
CABLE_forest_mask <- CABLE_forest_mask [, length (CABLE_forest_mask [1, ]):1]

# Read phenology mask
#----------------------------------------------------------------------------------------
CABLE_pheno_mask    <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/phenology/CABLE-POP_cruncep_lai_annual_1901_2015_phenology_mask.nc'
ncin <-  nc_open (nc_name)
CABLE_pheno_mask <- ncvar_get (ncin, 'phen_max_lai_phen_number')
CABLE_pheno_mask <- CABLE_pheno_mask [,length (CABLE_pheno_mask [1, ]):1]

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------
CABLE_CRUN_Cveg <- array (NA, dim = c (720, 360))
nc_name  <- 'data/CABLE-POP/CABLE-POP_cruncep_cveg_year_1901_2015.nc4'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'cveg', start = c (1, 1, 1, 85), count = c (720, 360, 10, 30)) # Extract variable
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the common period
CABLE_CRUN_Cveg <- tmp.array2 [,length (tmp.array2 [1, ]):1]

# Read CRUN NPP layer
#----------------------------------------------------------------------------------------
CABLE_CRUN_NPP <- array (NA, dim = c (720, 360))
nc_name  <- 'data/CABLE-POP/CABLE-POP_cruncep_npp_month_1901_2015.nc4'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'npp', start = c (1, 1, 1, 85*12+1), count = c (720, 360, 10, 30*12)) # Extract variable
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the 
CABLE_CRUN_NPP <- tmp.array2 [,length (tmp.array2 [1, ]):1] * 60.0 * 60.0 * 24.0 * 365.233333 # 7 /30 years in that period are leap years 

# Calculate CRUN tau_NPP for CABLE
#----------------------------------------------------------------------------------------
CABLE_CRUN_nta <- CABLE_CRUN_Cveg / CABLE_CRUN_NPP
CABLE_CRUN_nta [is.na (CABLE_forest_mask) | CABLE_forest_mask == 2] <- NA

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------
CABLE_CRUN_nta_adj <- CABLE_CRUN_Cveg / (CABLE_CRUN_NPP * (forest_mask_NA / 100.0))

# Clean up
#----------------------------------------------------------------------------------------
rm (tmp.array, tmp.array1, tmp.array2)
#========================================================================================