#========================================================================================
# Read LPJmL Cveg and mortality to calculate the actual turnover. The mortality
# components are all explained in README_totalmortflux.txt
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
nc_name <- 'forest_mask/hansen_forested_frac_05.nc4'
ncin <-  nc_open (nc_name)
forest_mask <- ncvar_get (ncin, 'forested_50_percent')
forest_mask_NA <- forest_mask 
forest_mask_NA [forest_mask_NA < threshold] <- NA

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

# Initialise CRUN LPJmL_mort = cmort_fire + cmort_growth_efficiency + 
#                              cmort_neg_allocation + cmort_neg_biomass + 
#                              cmort_heatstress + cmort_shade
#----------------------------------------------------------------------------------------
LPJmL_CRUN_mort <- array (NA, dim = c (720, 360))

# Get cmort_fire [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'lpjml_cruncep_cmort_fire_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp <- ncvar_get (ncin, 'cmort_fire', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract variable

# Get cmort_groweff [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'lpjml_cruncep_cmort_growth_efficiency_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp1 <- ncvar_get (ncin, 'cmort_growth_efficiency', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract variable

# Get cmort_neg_allocation [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'lpjml_cruncep_cmort_neg_allocation_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp2 <- ncvar_get (ncin, 'cmort_neg_allocation', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract variable

# Get cmort_neg_biomass [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'lpjml_cruncep_cmort_neg_biomass_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp3 <- ncvar_get (ncin, 'cmort_neg_biomass', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract variable

# Get cmort_heatstress [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'lpjml_cruncep_cmort_heatstress_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp4 <- ncvar_get (ncin, 'cmort_heatstress', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract variable

# Get cmort_shade [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'lpjml_cruncep_cmort_shade_annual_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp5 <- ncvar_get (ncin, 'cmort_shade', start = c (1, 1, 1, 85), count = c (720, 279, 9, 30)) # Extract variable

# Sum the mort components and delete them
#----------------------------------------------------------------------------------------
tmp.array <- tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5
rm (tmp, tmp1, tmp2, tmp3, tmp4, tmp5)

# Sum across PFTs
#----------------------------------------------------------------------------------------
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
rm (tmp.array)

# Average across the common period time
#----------------------------------------------------------------------------------------
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the common period
tmp.array3 <- cbind (matrix (data = NA, ncol = 67, nrow = 720), tmp.array2, matrix (data = NA, ncol = 14, nrow = 720)) # Add lines to get full North South extent
tmp.array3 <- tmp.array3 * 60.0 * 60.0 * 24.0 * 365.233333
rm (tmp.array1, tmp.array2)

# Clean-up
#----------------------------------------------------------------------------------------
LPJmL_CRUN_mort <- tmp.array3
rm (tmp.array3)

# Calculate tau
#----------------------------------------------------------------------------------------
LPJmL_CRUN_tau <- LPJmL_CRUN_Cveg / LPJmL_CRUN_mort
LPJmL_CRUN_tau [is.na (LPJmL_forest_mask) | LPJmL_forest_mask == 2] <- NA
#========================================================================================
