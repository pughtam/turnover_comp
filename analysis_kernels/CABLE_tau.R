#========================================================================================
# Read CABLE Cveg and mortality to calculate the actual turnover. The mortality
# components are all explained in README_totalmortflux.txt.
#----------------------------------------------------------------------------------------

# Load dependencies
#----------------------------------------------------------------------------------------
if (!exists ('nc_open',   mode = 'function')) library (ncdf4)

# Threshold for the forest cover
#----------------------------------------------------------------------------------------
threshold <- 10.0

# Read the forest mask
#----------------------------------------------------------------------------------------
forest_mask    <- array (NA, dim = c (720, 360))
forest_mask_NA <- array (NA, dim = c (720, 360))
nc_name <- 'hansen_forested_frac_05.nc4'
ncin <-  nc_open (nc_name)
forest_mask <- ncvar_get (ncin, 'forested_50_percent')
forest_mask_NA <- forest_mask 
forest_mask_NA [forest_mask_NA < threshold] <- NA

# Read forest mask
#----------------------------------------------------------------------------------------
CABLE_forest_mask <- array (NA, dim = c (720, 360))
nc_name <- 'CABLE-POP_cruncep_lai_annual_1901_2015_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
CABLE_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
CABLE_forest_mask [CABLE_forest_mask == 2] <- NA
CABLE_forest_mask <- CABLE_forest_mask [, length (CABLE_forest_mask [1, ]):1]

# Read phenology mask
#----------------------------------------------------------------------------------------
CABLE_pheno_mask <- array (NA, dim = c (720, 360))
nc_name <- 'CABLE-POP_cruncep_lai_annual_1901_2015_phenology_mask.nc'
ncin <-  nc_open (nc_name)
CABLE_pheno_mask <- ncvar_get (ncin, 'phen_max_lai_phen_number')
CABLE_pheno_mask <- CABLE_pheno_mask [,length (CABLE_pheno_mask [1, ]):1]

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------
CABLE_CRUN_Cveg <- array (NA, dim = c (720, 360))
nc_name  <- 'CABLE-POP_cruncep_cveg_year_1901_2015.nc4'
ncin <-  nc_open (nc_name)
tmp.array <- ncvar_get (ncin, 'cveg', start = c (1, 1, 1, 85), count = c (720, 360, 10, 30)) # Extract variable
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the 
CABLE_CRUN_Cveg <- tmp.array2 [,length (tmp.array2 [1, ]):1]

# Initialise CRUN CABLE_mort = cmortwoodresourcelim + 
#                              cmortwooddist + 
#                              cmortwoodresourcecrowding
#----------------------------------------------------------------------------------------
CABLE_CRUN_mort <- array (NA, dim = c (720, 360))

# Get cmortwoodresourcelim [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'CABLE-POP_cruncep_cmortwoodresourcelim_month_1901_2015.nc4'
ncin <-  nc_open (nc_name)
tmp <- ncvar_get (ncin, 'cmortwoodresourcelim', start = c (1, 1, 1, 85*12+1), count = c (720, 360, 10, 30*12)) # Extract variable

# Get cmortwooddist [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'CABLE-POP_cruncep_cmortwooddist_month_1901_2015.nc4'
ncin <-  nc_open (nc_name)
tmp1 <- ncvar_get (ncin, 'cmortwooddist', start = c (1, 1, 1, 85*12+1), count = c (720, 360, 10, 30*12)) # Extract variable

# Get cmortwoodresourcecrowding [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'CABLE-POP_cruncep_cmortwoodresourcecrowding_month_1901_2015.nc4'
ncin <-  nc_open (nc_name)
tmp2 <- ncvar_get (ncin, 'cmortwoodresourcecrowding', start = c (1, 1, 1, 85*12+1), count = c (720, 360, 10, 30*12)) # Extract variable

# Sum the mort components and delete them
#----------------------------------------------------------------------------------------
tmp.array <- tmp + tmp1 + tmp2

# Delete temporary variables
#----------------------------------------------------------------------------------------
rm (tmp1, tmp2, tmp)

# Sum across PFTs
#----------------------------------------------------------------------------------------
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
rm (tmp.array)

# Average across the common period time
#----------------------------------------------------------------------------------------
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the common period
tmp.array2 <- tmp.array2 [,length (tmp.array2 [1, ]):1] * 60.0 * 60.0 * 24.0 * 365.233333
rm (tmp.array1)

# Clean-up
#----------------------------------------------------------------------------------------
CABLE_CRUN_mort <- tmp.array2
rm (tmp.array2)

# Calculate tau
#----------------------------------------------------------------------------------------
CABLE_CRUN_tau <- CABLE_CRUN_Cveg / CABLE_CRUN_mort
CABLE_CRUN_tau [is.na (CABLE_forest_mask) | CABLE_forest_mask == 2] <- NA
#========================================================================================