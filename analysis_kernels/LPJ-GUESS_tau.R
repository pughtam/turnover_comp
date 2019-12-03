#========================================================================================
# Read LPJ-GUESS Cveg and mortality to calculate the actual turnover. The mortality
# components are all explained in README_totalmortflux.txt.
#----------------------------------------------------------------------------------------

# Read LPJ_GUESS forest mask
#----------------------------------------------------------------------------------------
LPJ_GUESS_forest_mask <- array (NA, dim = c (720, 360))
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

# Initialise CRUN LPJ_GUESS_mort = cmort_fire + cmort_groweff + cmort_dist + 
#                                  cmort_badallom + cmort_bioclim + 
#                                  cmort_mmin + cmort_nbiom
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_mort <- array (NA, dim = c (720, 360))

# Get cmort_fire [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/LPJ-GUESS/CRUNCEP2/lpj-guess_cruncep_cmort_fire_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp <- ncvar_get (ncin, 'cmort_fire', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract variable

# Get cmort_groweff [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/LPJ-GUESS/CRUNCEP2/lpj-guess_cruncep_cmort_groweff_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp1 <- ncvar_get (ncin, 'cmort_groweff', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract variable

# Get cmort_dist [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/LPJ-GUESS/CRUNCEP2/lpj-guess_cruncep_cmort_dist_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp2 <- ncvar_get (ncin, 'cmort_dist', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract variable

# Get cmort_badallom [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/LPJ-GUESS/CRUNCEP2/lpj-guess_cruncep_cmort_badallom_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp3 <- ncvar_get (ncin, 'cmort_badallom', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract variable

# Get cmort_bioclim [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/LPJ-GUESS/CRUNCEP2/lpj-guess_cruncep_cmort_bioclim_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp4 <- ncvar_get (ncin, 'cmort_bioclim', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract variable

# Get cmort_mmin [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/LPJ-GUESS/CRUNCEP2/lpj-guess_cruncep_cmort_mmin_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp5 <- ncvar_get (ncin, 'cmort_mmin', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract variable

# Get cmort_nbiom [kg m-2 s-1]
#----------------------------------------------------------------------------------------
nc_name  <- 'data/LPJ-GUESS/CRUNCEP2/lpj-guess_cruncep_cmort_nbiom_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
tmp6 <- ncvar_get (ncin, 'cmort_nbiom', start = c (1, 1, 1, 85), count = c (720, 360, 11, 30)) # Extract variable

# Sum the mort components and delete them
#----------------------------------------------------------------------------------------
tmp.array <- tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6
rm (tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

# Sum across PFTs
#----------------------------------------------------------------------------------------
tmp.array1 <- apply (tmp.array, c (1, 2, 4), sum, na.rm = T) # Sum variable across pfts
rm (tmp.array)

# Average across the common period time
#----------------------------------------------------------------------------------------
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the common period
tmp.array2 <- tmp.array2 * 60.0 * 60.0 * 24.0 * 365.233333
rm (tmp.array1)

# Clean-up
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_mort <- tmp.array2
rm (tmp.array2)

# Calculate tau
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_tau <- LPJ_GUESS_CRUN_Cveg / LPJ_GUESS_CRUN_mort
LPJ_GUESS_CRUN_tau [is.na (LPJ_GUESS_forest_mask) | LPJ_GUESS_forest_mask == 2] <- NA

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------
LPJ_GUESS_CRUN_tau_adj <- LPJ_GUESS_CRUN_tau * (forest_mask_NA / 100.0)

#========================================================================================
