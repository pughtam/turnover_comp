#==============================================================================#
# Read SEIB Cveg and mortality to calculate the actual turnover. The mortality
# components are all explained in README_totalmortflux.txt
#------------------------------------------------------------------------------#

# Load dependencies
#------------------------------------------------------------------------------#
if (!exists ('nc_open',   mode = 'function')) library (ncdf4)
if (!exists ('pl.global', mode = 'function')) source ('R_scripts//image_functions.R')

# Threshold for the forest cover
#----------------------------------------------------------------------------------------#
threshold <- 10.0

# Read the forest mask
#----------------------------------------------------------------------------------------#
forest_mask    <- array (NA, dim = c (720, 360))
forest_mask_NA <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/forest_mask/hansen_forested_frac_05.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
forest_mask <- ncvar_get (ncin, 'forested_50_percent')
forest_mask_NA <- forest_mask 
forest_mask_NA [forest_mask_NA < threshold] <- NA
#pl.global (forest_mask_NA, ulimit = 100)

# Read forest mask
#----------------------------------------------------------------------------------------#
SEIB_forest_mask    <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/forest_mask/seib_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'
ncin <-  nc_open (nc_name)
#print (ncin)
SEIB_forest_mask <- ncvar_get (ncin, 'forest_30yr_any_10_years')
SEIB_forest_mask [SEIB_forest_mask == 2] <- NA
SEIB_forest_mask <- SEIB_forest_mask [,length (SEIB_forest_mask [1, ]):1]
#pl.global (SEIB_forest_mask, ulimit = 1)

# Read phenology mask
#----------------------------------------------------------------------------------------#
SEIB_pheno_mask  <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/phenology/seib_cruncep_lai_annual_1901_2014_phenology_mask.nc'
ncin <-  nc_open (nc_name)
#print (ncin)
SEIB_pheno_mask <- ncvar_get (ncin, 'phen_max_lai_phen_number')
SEIB_pheno_mask <- SEIB_pheno_mask [, length (SEIB_pheno_mask [1, ]):1]
#pl.global (SEIB_pheno_mask)

# Read CRUN Cveg layer
#----------------------------------------------------------------------------------------#
SEIB_CRUN_Cveg <- array (NA, dim = c (720, 360))
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cveg_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp.array <- ncvar_get (ncin, 'cveg', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract only the years 1985 to 2014
tmp.array1 <- apply (tmp.array, c (1, 2, 3), sum, na.rm = T) # Sum variable across pfts
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over time 
tmp.array3 <- tmp.array2 [, length (tmp.array2 [1, ]):1] # Turn the world up-side-down
tmp.array3 [is.na (forest_mask_NA) | is.na (SEIB_forest_mask)] <- NA # Apply the forest mask
SEIB_CRUN_Cveg <- tmp.array3
#pl.global (SEIB_CRUN_Cveg)

# Initialise CRUN SEIB_mort = cmort_fire_atmosp + cmort_fire_litter + 
#                             cmort_greff + cmort_gap + cmort_heat + 
#                             cmort_lim + cmort_etc
SEIB_CRUN_mort <- array (NA, dim = c (720, 360))

# Get cmort_fire_atmosp [kg m-2 s-1]
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cmort_fire_atmosp_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp <- ncvar_get (ncin, 'cmort_fire_atmosp', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract variable

# Get cmort_fire_litter [kg m-2 s-1]
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cmort_fire_litter_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp1 <- ncvar_get (ncin, 'cmort_fire_litter', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract variable

# Get cmort_greff [kg m-2 s-1]
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cmort_greff_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp2 <- ncvar_get (ncin, 'cmort_greff', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract variable

# Get cmort_gap [kg m-2 s-1]
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cmort_gap_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp3 <- ncvar_get (ncin, 'cmort_gap', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract variable

# Get cmort_heat [kg m-2 s-1]
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cmort_heat_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp4 <- ncvar_get (ncin, 'cmort_heat', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract variable

# Get cmort_lim [kg m-2 s-1]
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cmort_lim_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp5 <- ncvar_get (ncin, 'cmort_lim', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract variable

# Get cmort_etc [kg m-2 s-1]
nc_name  <- 'data/SEIB/cruncep_netcdfs_v6/seib_cruncep_cmort_etc_annual_1901_2014.nc4'
ncin <-  nc_open (nc_name)
#print (ncin)
tmp6 <- ncvar_get (ncin, 'cmort_etc', start = c (1, 1, 85, 1), count = c (720, 360, 30, 14)) # Extract variable

# Sum all components
tmp.array <- tmp + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6
rm (tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

# Sum across PFTs
tmp.array1 <- apply (tmp.array, c (1, 2, 3), sum, na.rm = T) # Sum variable across pfts
rm (tmp.array)

# Average across the common period time
tmp.array2 <- apply (tmp.array1, c (1, 2), mean, na.rm = T) # Calculate mean over the common period
tmp.array3 <- tmp.array2 [, length (tmp.array2 [1, ]):1]
tmp.array4 <- tmp.array3 * 60.0 * 60.0 * 24.0 * 365.233333 
rm (tmp.array1, tmp.array2, tmp.array3)

# Clean-up
SEIB_CRUN_mort <- tmp.array4
#pl.global (SEIB_CRUN_mort)
rm (tmp.array4)

# Calculate tau
SEIB_CRUN_tau <- SEIB_CRUN_Cveg / SEIB_CRUN_mort
SEIB_CRUN_tau [is.na (SEIB_forest_mask) | SEIB_forest_mask == 2] <- NA
SEIB_CRUN_tau1 <- SEIB_CRUN_tau
SEIB_CRUN_tau2 <- SEIB_CRUN_tau
SEIB_CRUN_tau1 [SEIB_CRUN_tau > 200] <- 200
SEIB_CRUN_tau2 [SEIB_CRUN_tau > 200] <- NA
#pl.global (SEIB_CRUN_tau, ulimit = 100)

# Weight residence time by forest fraction
#----------------------------------------------------------------------------------------#
SEIB_CRUN_tau_adj <- SEIB_CRUN_tau * (forest_mask_NA / 100.0)
SEIB_CRUN_tau1_adj <- SEIB_CRUN_tau1 * (forest_mask_NA / 100.0)
SEIB_CRUN_tau2_adj <- SEIB_CRUN_tau2 * (forest_mask_NA / 100.0)
#pl.global (SEIB_CRUN_tau_adj, ulimit = 100)

# plot density kernel
#----------------------------------------------------------------------------------------#
d_SEIB_nta   <- density (SEIB_CRUN_nta_adj,  na.rm = T)
d_SEIB_tau   <- density (SEIB_CRUN_tau_adj,  na.rm = T)
d_SEIB_tau1  <- density (SEIB_CRUN_tau1_adj, na.rm = T)
d_SEIB_tau2  <- density (SEIB_CRUN_tau2_adj, na.rm = T)

if (PLOT){
  colours <- c ('#99999999', '#106470', '#B60026', '#55A51C', '#EA7125', '#8F2BBC', '#91b9a4', '#68ACE5')
  
  par (mfrow = c (1, 1))
  par (mar = c (5, 5, 1, 1))
  plot (d_SEIB_tau,
        xlab = 'turnover time (years)',
        ylab = 'density',
        col = colours [8],
        lwd = 2,
        main = '',
        xlim = c (0, 85),
        ylim = c (-0.01, 0.07))
  lines (d_SEIB_tau1, 
         lwd = 2, 
         lty = 2,
         col = colours [8])
  lines (d_SEIB_tau2, 
         lwd = 2, 
         lty = 3,
         col = colours [8])
  abline (h = 0, lwd = 1, col = 'black')
  lines (d_SEIB_nta, 
         lwd = 2, 
         lty = 2,
         col = colours [8])
  text (x = 80,
        y = 0.07,
        labels = 'SEIB',
        pos = 1)
  legend (x = 50,
          y = 0.05,
          col = colours [8],
          lty = c (1, 2, 3),
          lwd = 2,
          bg = 'transparent',
          box.lty = 0,
          legend = c ('with outliers','outliers = 200', 'outliers = NA'))
          #legend = c (expression (tau), expression (tau[eq])))
}
#==============================================================================#
