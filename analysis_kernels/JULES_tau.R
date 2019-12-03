#========================================================================================
# Read JULES Cveg and mortality to calculate the actual turnover. The mortality
# components are all explained in README_totalmortflux.txt.
#----------------------------------------------------------------------------------------

# Read JULES forest mask
#----------------------------------------------------------------------------------------
JULES_forest_mask <- array (NA, dim = c (192, 112+29+4))
nc_name <- 'Shared_files/forest_mask/JULESC2_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
tmp.array1 <- ncvar_get (ncin, 'forest_30yr_any_10_years')
tmp.array2 <- cbind (matrix (rep (NA, 192 * 29), nrow = 192), rbind (tmp.array1 [97:192, ], tmp.array1 [1:96, ]), matrix (rep (NA, 192 * 4), nrow = 192))
JULES_forest_mask <- tmp.array2
JULES_forest_mask [JULES_forest_mask == 2] <- NA

# Read JULES phenology mask
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

# Initialise CRUN JULES_mort = wood_litC + lit_c_dyn_wood + lit_c_dyn_leaf + 
#                              lit_c_dyn_root
#----------------------------------------------------------------------------------------
JULES_CRUN_mort <- array (NA, dim = c (192, 112+29+4))

# Get wood_litC [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/JULES/cruncep2/JULESC2_CRUNCEP_wood_litC_Monthly_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp <- ncvar_get (ncin, 'wood_litC', start =  c (1, 1, 1, 85*12+1), count = c (192, 112, 9, 12*29)) # Extract variable

# Get lit_C_dyn_wood [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/JULES/cruncep2/JULESC2_CRUNCEP_lit_c_dyn_wood_Monthly_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp1 <- ncvar_get (ncin, 'lit_c_dyn_wood', start =  c (1, 1, 1, 85*12+1), count = c (192, 112, 9, 12*29)) # Extract variable

# Get lit_C_dyn_leaf [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/JULES/cruncep2/JULESC2_CRUNCEP_lit_c_dyn_leaf_Monthly_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp2 <- ncvar_get (ncin, 'lit_c_dyn_leaf', start =  c (1, 1, 1, 85*12+1), count = c (192, 112, 9, 12*29)) # Extract variable

# Get lit_C_dyn_root [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/JULES/cruncep2/JULESC2_CRUNCEP_lit_c_dyn_root_Monthly_1901_2014.nc'
ncin <-  nc_open (nc_name)
tmp3 <- ncvar_get (ncin, 'lit_c_dyn_root', start =  c (1, 1, 1, 85*12+1), count = c (192, 112, 9, 12*29)) # Extract variable

# Sum the mort components and delete them
#----------------------------------------------------------------------------------------
tmp.array <- tmp + tmp1 + tmp2 + tmp3
rm (tmp, tmp1, tmp2, tmp3)

# Sum across PFTs
#----------------------------------------------------------------------------------------
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
rm (tmp.array)

# Average across the common period time
#----------------------------------------------------------------------------------------
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the common period
tmp.array2 <- cbind (matrix (rep (NA, 192 * 29), nrow = 192), rbind (tmp.array2 [97:192, ], tmp.array2 [1:96, ]), matrix (rep (NA, 192 * 4), nrow = 192))
tmp.array2 <- tmp.array2 * 60.0 * 60.0 * 24.0 * 365.233333
rm (tmp.array1)

# Clean-up
#----------------------------------------------------------------------------------------
JULES_CRUN_mort <- tmp.array2
rm (tmp.array2)

# Calculate tau
#----------------------------------------------------------------------------------------
JULES_CRUN_tau <- JULES_CRUN_Cveg / JULES_CRUN_mort
JULES_CRUN_tau [is.na (JULES_forest_mask) | JULES_forest_mask == 2] <- NA

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------
JULES_CRUN_tau_adj <- JULES_CRUN_tau * (forest_mask_julesgrid_NA / 100.0)

#========================================================================================
