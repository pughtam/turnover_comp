#========================================================================================
# This script calcualtes the forest turnover times for ORCHIDEE.
#----------------------------------------------------------------------------------------

# Read ORCHIDEE forest mask
#----------------------------------------------------------------------------------------
ORCHIDEE_forest_mask    <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/forest_mask/ORCHIDEE_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
ORCHIDEE_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
ORCHIDEE_forest_mask [ORCHIDEE_forest_mask == 2] <- NA
ORCHIDEE_forest_mask <- ORCHIDEE_forest_mask [, length (ORCHIDEE_forest_mask [1, ]):1]

# Read phenology mask
#----------------------------------------------------------------------------------------
ORCHIDEE_pheno_mask  <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/phenology/ORCHIDEE_cruncep_lai_annual_1901_2014_phenology_mask.nc'
ncin <-  nc_open (nc_name)
ORCHIDEE_pheno_mask <- ncvar_get (ncin, 'phen_max_lai_phen_number')
ORCHIDEE_pheno_mask <- ORCHIDEE_pheno_mask [, length (ORCHIDEE_pheno_mask [1, ]):1]

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_Cveg <- array (NA, dim = c (720, 360))
nc_name  <- 'data/ORCHIDEE/cruncep2/ORCHIDEEr3085_cruncep_cVeg_13pft_year_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'cVeg', star = c (1, 1, 1, 85), count = c (720, 360, 13, 30)) # Extract only the years 1982 to 2011
tmp.array [tmp.array >= 1.0e+6] <- NA
tmp.array1 <- apply (tmp.array, c (1, 2, 3), mean, na.rm = T) # Calculate mean over the 

# Read vegetation fraction
#----------------------------------------------------------------------------------------
nc_name  <- 'data/ORCHIDEE/cruncep2/ORCHIDEEr3085_cruncep_landCoverFrac_13pft_year_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp <- ncvar_get (ncin, 'landCoverFrac', star = c (1, 1, 1, 85), count = c (720, 360, 13, 30)) # Extract only the years 1982 to 2011
tmp [tmp >= 1.0e+6] <- NA
tmp1 <- apply (tmp, c (1, 2, 3), mean, na.rm = T) # Calculate mean over time

# Read continental fraction
#----------------------------------------------------------------------------------------
nc_name  <- 'data/ORCHIDEE/cruncep2/ORCHIDEEr3085_cruncep_contfrac_gridcell_year_1901_2014.nc'
ncin <-  nc_open (nc_name)
contfrac <- ncvar_get (ncin, 'contfrac', star = c (1, 1, 85), count = c (720, 360, 30)) # Extract only the years 1982 to 2011
contfrac [contfrac >= 1.0e+6] <- NA
contfrac <- apply (contfrac, c (1, 2), mean, na.rm = T) # Calculate mean over the 

# Adjust each pfts cVeg by the land cover fraction
#----------------------------------------------------------------------------------------
tmp.array2 <- tmp.array1
for (i in 1:13) {
  tmp.array2 [,, i] <- (tmp1 [,, i] / 100.0) * tmp.array1 [, , i] * contfrac
}
tmp.array3 <- apply (tmp.array2, c (1, 2), sum, na.rm = T) # Sum variable across pfts
tmp.array4 <- tmp.array3 [, length (tmp.array3 [1, ]):1] # Turn the world up-side-down
tmp.array4 [is.na (forest_mask_NA) | is.na (ORCHIDEE_forest_mask)] <- NA # Apply the forest mask

ORCHIDEE_CRUN_Cveg <- tmp.array4

# Read CRUN NPP layer
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_NPP <- array (NA, dim = c (720, 360))
nc_name  <- 'data/ORCHIDEE/cruncep2/ORCHIDEEr3085_cruncep_npp_13pft_year_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'npp', star = c (1, 1, 1, 85), count = c (720, 360, 13, 30))# Extract only the years 1982 to 2011
tmp.array [tmp.array >= 1.0e+6] <- NA

# Read veget_max
#----------------------------------------------------------------------------------------
nc_name  <- 'data/ORCHIDEE/cruncep2/ORCHIDEEr3085_cruncep_VEGET_MAX_year_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp <- ncvar_get (ncin, 'VEGET_MAX', star = c (1, 1, 1, 85), count = c (720, 360, 13, 30))# Extract only the years 1982 to 2011
tmp [tmp >= 1.0e+6] <- NA
tmp.array <- tmp * tmp.array

tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the 
tmp.array3 <- tmp.array2 [, length (tmp.array2 [1, ]):1] # Turn the world up-side-down
tmp.array3 [is.na (forest_mask_NA) | is.na (ORCHIDEE_forest_mask)] <- NA
ORCHIDEE_CRUN_NPP <- tmp.array3 * 60 * 60 * 24 * 365.2333333

# Calculate CRUN tau_NPP for ORCHIDEE
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_nta <- ORCHIDEE_CRUN_Cveg / ORCHIDEE_CRUN_NPP
ORCHIDEE_CRUN_nta [is.na (ORCHIDEE_forest_mask) | ORCHIDEE_forest_mask == 2] <- NA

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_nta_adj <- ORCHIDEE_CRUN_nta * (forest_mask_NA / 100.0)

# Clean up
#----------------------------------------------------------------------------------------#
rm (tmp.array, tmp.array1, tmp.array2)
#========================================================================================