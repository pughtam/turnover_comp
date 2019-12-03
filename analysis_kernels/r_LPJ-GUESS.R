#========================================================================================
# This script calcualtes the forest turnover times for LPJ-GUESS.
#----------------------------------------------------------------------------------------

# Read LPJ-GUESS forest mask
#----------------------------------------------------------------------------------------
LPJ_GUESS_forest_mask    <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files//forest_mask/lpj-guess_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
LPJ_GUESS_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
LPJ_GUESS_forest_mask [LPJ_GUESS_forest_mask == 2] <- NA

# Read phenology mask
#----------------------------------------------------------------------------------------
LPJ_GUESS_pheno_mask  <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/phenology/lpj-guess_cruncep_lai_annual_1901_2014_phenology_mask.nc'
ncin <-  nc_open (nc_name)
LPJ_GUESS_pheno_mask <- ncvar_get (ncin, 'phen_max_lai_phen_number')

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_Cveg <- array (NA, dim = c (720, 360))
nc_name  <- 'data/LPJ-GUESS/CRUNCEP/lpj-guess_cruncep_cveg_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'cveg', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract only the years 1982 to 2011
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the 
tmp.array2 [is.na (LPJ_GUESS_forest_mask)] <- NA
LPJ_GUESS_CRUN_Cveg <- tmp.array2

# Read CRUN NPP layer
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_NPP <- array (NA, dim = c (720, 360))
nc_name  <- 'data/LPJ-GUESS/CRUNCEP/lpj-guess_cruncep_npp_monthly_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'npp', start = c (1, 1, 1, 12*84+1), count = c (720, 360, 11, 12*30)) # Extract only the years 1982 to 2011
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the 
tmp.array2 [is.na (forest_mask_NA) | is.na (LPJ_GUESS_forest_mask)] <- NA
LPJ_GUESS_CRUN_NPP <- tmp.array2 * 60.0 * 60.0 * 24.0 * 365.233333 # 7 /30 years in that period are leap years 

# Calculate CRUN LPJ-GUESS tau_NPP
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_nta <- LPJ_GUESS_CRUN_Cveg / LPJ_GUESS_CRUN_NPP
LPJ_GUESS_CRUN_nta [is.na (LPJ_GUESS_forest_mask) | LPJ_GUESS_forest_mask == 2] <- NA

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_nta_adj <- LPJ_GUESS_CRUN_nta * (forest_mask_NA / 100.0)

# Clean up
#----------------------------------------------------------------------------------------
rm (tmp.array, tmp.array1, tmp.array2)
#========================================================================================