#========================================================================================
# Read ORCHIDEE Cveg and mortality to calculate the actual turnover. The mortality
# components are all explained in README_totalmortflux.txt.
#----------------------------------------------------------------------------------------

# Read ORCHIDEE forest mask
#----------------------------------------------------------------------------------------
ORCHIDEE_forest_mask <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/forest_mask/ORCHIDEE_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
ORCHIDEE_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
ORCHIDEE_forest_mask [ORCHIDEE_forest_mask == 2] <- NA
ORCHIDEE_forest_mask <- ORCHIDEE_forest_mask [, length (ORCHIDEE_forest_mask [1, ]):1]

# Read phenology mask
#----------------------------------------------------------------------------------------
ORCHIDEE_pheno_mask <- array (NA, dim = c (720, 360))
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

# Initialise CRUN ORCHIDEE_mort = total_turn
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_mort <- array (NA, dim = c (720, 360))

# Get total_turn [g m-2 d-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/ORCHIDEE/cruncep2/ORCHIDEEr3085_cruncep_total_turn_13pft_year_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp <- ncvar_get (ncin, 'total_turn', start = c (1, 1, 1, 85), count = c (720, 360, 13, 30)) # Extract variable
tmp [tmp >= 1.0e+6] <- NA
tmp.array <- apply (tmp, c (1, 2, 3), mean, na.rm = T) # Calculate mean over the common period

# Multiply by vegfrac and contfrac
#----------------------------------------------------------------------------------------
for (i in 1:13) {
  tmp.array1 [,, i] <- (tmp1 [,, i] / 100.0) * tmp.array [, , i] * contfrac
}

# Sum across PFTs
#----------------------------------------------------------------------------------------
tmp.array2 <- apply (tmp.array1, c (1, 2), sum, na.rm = T) # Sum variable across pfts
rm (tmp)

# Average across the common period time
#----------------------------------------------------------------------------------------
tmp.array3 <- tmp.array2 [, length (tmp.array2 [1, ]):1]
tmp.array4 <- tmp.array3 * 365.233333 / 1000.0
rm (tmp.array1, tmp.array2, tmp.array3)

# Clean-up
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_mort <- tmp.array4
rm (tmp.array4)

# Calculate tau
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_tau <- ORCHIDEE_CRUN_Cveg / ORCHIDEE_CRUN_mort
ORCHIDEE_CRUN_tau [is.na (ORCHIDEE_forest_mask) | ORCHIDEE_forest_mask == 2] <- NA

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------
ORCHIDEE_CRUN_tau_adj <- ORCHIDEE_CRUN_tau * (forest_mask_NA / 100.0)

#========================================================================================