#========================================================================================
# Plot figures 3 and 6 for Pugh et al. (2019) titled 
# "Understanding the uncertainty in global forest carbon turnover".
#
# This code was writen by Tim T. Rademacher (trademacher@fas.harvard.edu).
#
# Last modified: 2019/12/03
#----------------------------------------------------------------------------------------

# Load dependencies
#----------------------------------------------------------------------------------------
library ('R.matlab') 
library ('ncdf4')

# Threshold for the definition of forest in grid cell forest cover (%)
#----------------------------------------------------------------------------------------
threshold <- 10.0

# Read the forest mask
#----------------------------------------------------------------------------------------
forest_mask    <- array (NA, dim = c (720, 360))
forest_mask_NA <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/forest_mask/hansen_forested_frac_05.nc4'
ncin <-  nc_open (nc_name)
forest_mask <- ncvar_get (ncin, 'forested_50_percent')
forest_mask_NA <- forest_mask 
forest_mask_NA [forest_mask_NA < threshold] <- NA

# Read the canopy cover mask
#----------------------------------------------------------------------------------------
canopy_mask    <- array (NA, dim = c (720, 360))
canopy_mask_NA <- array (NA, dim = c (720, 360))
nc_name <- 'Shared_files/forest_mask/hansen_canopy_cover_frac.nc4'
ncin <-  nc_open (nc_name)
canopy_mask <- ncvar_get (ncin, 'canopy_cover')
canopy_mask_NA <- canopy_mask 
canopy_mask_NA [canopy_mask_NA < threshold] <- NA

# Read the canopy cover mask for JULES (not 0.5 degree x 0.5 degree resolution)
#----------------------------------------------------------------------------------------
nc_name <- 'Shared_files/forest_mask/hansen_forested_frac_julesgrid.nc4'
ncin <-  nc_open (nc_name)
forest_mask_julesgrid <- ncvar_get (ncin, 'forested_50_percent')
forest_mask_julesgrid_NA <- forest_mask_julesgrid [, length (forest_mask_julesgrid [1, ]):1]
forest_mask_julesgrid_NA [forest_mask_julesgrid_NA < threshold] <- NA

# Read the phenology mask
#----------------------------------------------------------------------------------------
nc_name <- 'Shared_files/phenology/ESA_forest_9regions_v2.nc'
ncin <-  nc_open (nc_name)
pheno_mask <- ncvar_get (ncin, 'region_mask')

# Read Nuno's vegetation carbon layer
#----------------------------------------------------------------------------------------
Cveg <- check.for.mod.var ('C_v1', 'carvalhais14')
Cveg [is.na (forest_mask_NA)] <- NA

# Read MODIS NPP layer
#----------------------------------------------------------------------------------------
NPP <- check.for.mod.var ('NPP', 'carvalhais14')
NPP [is.na (forest_mask_NA)] <- NA

# Calculate tau and weighted alternatives
#----------------------------------------------------------------------------------------
nta <- Cveg / NPP
nta_adj1 <- Cveg / (NPP * (forest_mask_NA / 100.0))
nta_adj2 <- Cveg / (NPP * (canopy_mask_NA / 100.0))

# Clean up to free up memory
#----------------------------------------------------------------------------------------
rm (nc_name, ncin, threshold)

# Set colour scheme for phenologies
#----------------------------------------------------------------------------------------#
colours <- c ('#91b9a4', '#106470', '#EB99A9', '#55A51C', '#EA7125', '#8F2BBC', '#D6083B')

# phenologies
#----------------------------------------------------------------------------------------#
# 1 = needleleaved evergreen
# 2 = needleleaved deciduous
# 3 = boreal broadleaved deciduous
# 4 = temperate broadleaved evergreen
# 5 = temperate broadleaved deciduous
# 6 = tropical broadleaved evergreen
# 7 = tropical broadleaved raingreen

# Adjust categories for observational phenology mask
#----------------------------------------------------------------------------------------#
pheno_mask_obs <- matrix (as.numeric (rep (NA, 720 * 360)), nrow = 720)
pheno_mask_obs [pheno_mask == 6] <- 1 # NE
pheno_mask_obs [pheno_mask == 7] <- 2 # ND
pheno_mask_obs [pheno_mask == 8] <- 3 # boreal BD or set it to 3
pheno_mask_obs [pheno_mask == 4] <- 4 # temperate BE
pheno_mask_obs [pheno_mask == 5] <- 5 # temperate BD
pheno_mask_obs [pheno_mask == 1] <- 6 # tropical BE
pheno_mask_obs [pheno_mask == 2] <- 7 # tropical BD
pheno_mask_obs [pheno_mask == 9] <- NA # tropical BD

# Function to add polygon to density kernel
#----------------------------------------------------------------------------------------
add.polygon <- function (x1, y1, y2, i = 1) {
  polygon (x   = c (x1, rev (x1)),
           y   = c (y1, rev (y2)),
           lty = 0,
           col = colours [i])
  
}
add.mean <- function (i, nta, pheno_mask) {
  points (x = mean (nta [pheno_mask == i], na.rm = T),
          y = -0.01,
          lwd = 2, pch = 19,
          col = colours [i])
}

# Function to plot density kernel by phenology for a dataset (e.g. TBM or satellite observations)
#----------------------------------------------------------------------------------------
pl.nta.pheno <- function (nta_in, pheno_mask_in, string, eq = T) {
  for (k in 1:7) {
    if (sum (nta_in [pheno_mask_in == k], na.rm = T) > 10) {
      assign (x = paste ('d_',as.character (k), sep = ''),
              value = density (nta_in [pheno_mask_in == k], na.rm = T, 
                               from = 0, to = 100))
    } else {
      assign (x = paste ('d_',as.character (k), sep = ''),
              value = NA)
    }
  }
  
  par (mfrow = c (1, 1))
  par (mar = c (5, 5, 1, 1))
  plot (d_1,
        xlab = expression (paste (tau[NPP],' (years)')),
        ylab = 'density',
        col = colours [1],
        lwd = 0,
        main = '',
        xlim = if (eq) {c (0, 35)} else  {c (0, 70)},
        ylim = if (eq) {c (-0.01, 0.16)} else {c (-0.01, 0.10)})
  abline (h = 0, lwd = 1, col = 'black')
  
  x1 <- d_1$x
  y1 <- rep (0, length (d_1$x))
  y2 <- rep (0, length (d_1$x))
  for (k in 1:7) {
    if (is.na (get (paste ('d_',as.character (k), sep = ''))) [1]) next
    y1 <- y2
    multiplier <- sum (pheno_mask_in == k, na.rm = T) / sum (!is.na (pheno_mask_in), na.rm = T)
    additor <- get (paste ('d_',as.character (k), sep = ''))$y
    y2 <- y2 + additor  * multiplier
    add.polygon (x1 = x1, y1 = y1, y2 = y2, i = k)
  }
  
  lines (x = x1,
         y = y2)
  abline (h = 0, lwd = 1, col = 'black')
  abline (h = -0.01, lwd = 1, col = 'grey')
  
  res <- sapply (1:7, add.mean, nta = nta_in, pheno_mask = pheno_mask_in); rm (res)
  
  text (x = 0,
        y = if (eq) {0.15} else {0.09},
        labels = string,
        pos = 4)
  legend (x = if (eq) {28} else {55},
          y = if (eq) {0.16} else {0.10},
          legend = c ('NE','ND','boreal BD','temperate BE','temperate BD',
                      'tropical BE','tropical BR'),
          col = colours,
          pch = 19,
          cex = 0.8,
          box.lty = 0)
}

# Read CABLE fields 
#----------------------------------------------------------------------------------------
source ('R_scripts/r_CABLE-POP.R')
rm (forest_mask_julesgrid, forest_mask_julesgrid_NA, forest_mask, forest_mask_NA, 
    GPP, NPP, Ceco, Csoil, Cveg)

# Plot density kernel for CABLE (panel (a) in Fig. 3)
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = CABLE_CRUN_nta_adj,
              pheno_mask = CABLE_pheno_mask,
              string = 'CABLE')

# Read JULES fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/r_JULES.R')
rm (forest_mask_julesgrid, forest_mask_julesgrid_NA, forest_mask, forest_mask_NA, 
    GPP, NPP, Ceco, Csoil, Cveg)

# Plot density kernel for JULES (panel (b) in Fig. 3)
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = JULES_CRUN_nta_adj,
              pheno_mask = JULES_pheno_mask,
              string = 'JULES')

# Read LPJ-GUESS fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/r_LPJ-GUESS.R')
rm (forest_mask_LPJ_GUESS, forest_mask_LPJ_GUESS_NA, forest_mask, forest_mask_NA, 
    GPP, NPP, Ceco, Csoil, Cveg)

# Plot density kernels for LPJ-GUESS (panel (c) in Fig. 3)
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = LPJ_GUESS_CRUN_nta_adj,
              pheno_mask = LPJ_GUESS_pheno_mask,
              string = 'LPJ-GUESS')

# Read LPJmL fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/r_LPJmL.R')
rm (forest_mask_LPJmL, forest_mask_LPJmL_NA, forest_mask, forest_mask_NA, 
    GPP, NPP, Ceco, Csoil, Cveg)

# Plot density kernels for LPJmL (panel (d) in Fig. 3)
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = LPJmL_CRUN_nta_adj,
              pheno_mask = LPJmL_pheno_mask,
              string = 'LPJmL')

# Read ORCHIDEE fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/r_ORCHIDEE.R')
rm (forest_mask_ORCHIDEE, forest_mask_ORCHIDEE_NA, forest_mask, forest_mask_NA, 
    GPP, NPP, Ceco, Csoil, Cveg)

# Plot density kernels for ORCHIDEE (panel (e) in Fig. 3)
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = ORCHIDEE_CRUN_nta_adj,
              pheno_mask = ORCHIDEE_pheno_mask,
              string = 'ORCHIDEE')

# Read SEIB fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/r_SEIB.R')
rm (forest_mask_SEIB, forest_mask_SEIB_NA, forest_mask, forest_mask_NA, 
    GPP, NPP, Ceco, Csoil, Cveg)

# Plot density kernels for SEIB (panel (f) in Fig. 3)
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = SEIB_CRUN_nta_adj,
              pheno_mask = SEIB_pheno_mask,
              string = 'SEIB')

# Plot observational fields (panel (g) in Fig. 3)
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = nta,
              pheno_mask = pheno_mask_obs,
              string = 'satellite-based')

# Plot figure 6 - mortality based turnover time (tau_mort / tau)
#----------------------------------------------------------------------------------------#

# Read CABLE fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/CABLE_tau.R')

# Plot density kernels for CABLE
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = CABLE_CRUN_tau_adj,
              pheno_mask = CABLE_pheno_mask,
              string = 'CABLE',
              eq = F)

# Read JULES fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/JULES_tau.R')

# Plot density kernels for JULES
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = JULES_CRUN_tau_adj,
              pheno_mask = JULES_pheno_mask,
              string = 'JULES',
              eq = F)

# Read LPJ-GUESS fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/LPJ-GUESS_tau.R')

# Plot density kernels for LPJ-GUESS
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = LPJ_GUESS_CRUN_tau_adj,
              pheno_mask = LPJ_GUESS_pheno_mask,
              string = 'LPJ-GUESS',
              eq = F)

# Read LPJmL fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/LPJmL_tau.R')

# Plot density kernels for LPJmL
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = LPJmL_CRUN_tau_adj,
              pheno_mask = LPJmL_pheno_mask,
              string = 'LPJmL',
              eq = F)

# Read ORCHIDEE fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/ORCHIDEE_tau.R')

# Plot density kernels for ORCHIDEE
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = ORCHIDEE_CRUN_tau_adj,
              pheno_mask = ORCHIDEE_pheno_mask,
              string = 'ORCHIDEE',
              eq = F)

# Read SEIB fields
#----------------------------------------------------------------------------------------#
source ('R_scripts/SEIB_tau.R')

# Plot density kernels for SEIB
#----------------------------------------------------------------------------------------#
pl.nta.pheno (nta = SEIB_CRUN_tau_adj,
              pheno_mask = SEIB_pheno_mask,
              string = 'SEIB',
              eq = F)

#========================================================================================#