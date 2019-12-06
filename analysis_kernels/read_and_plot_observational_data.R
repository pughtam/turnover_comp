#========================================================================================
# This script reads Nuno Carvalhais's vegetation biomass and Martin Thurner's NPP data to
# calcualte forest turnover times.
#----------------------------------------------------------------------------------------

# Require netcdf library
#----------------------------------------------------------------------------------------
if (!exists ('nc_open',   mode = 'function')) library (ncdf4)

# Threshold for the forest cover
#----------------------------------------------------------------------------------------
threshold <- 10.0

# Read the phenology mask
#----------------------------------------------------------------------------------------
nc_name <- 'ESA_forest_9regions.nc'
ncin <-  nc_open (nc_name)
pheno_mask <- ncvar_get (ncin, 'region_mask')
pheno_mask <- pheno_mask [, length (pheno_mask [1, ]):1]

# Read soil carbon and ecosystem carbon later layer from Carvalhais et al. (2014) to 
# calculate vegetation carbon layer.
#----------------------------------------------------------------------------------------
C_e <- read.table (file = 'C_e.txt', header = F)
C_e [C_e == -99] <- NA
C_s <- read.table (file = 'C_s.txt', header = F)
C_s [C_s == -99] <- NA
Cveg <- C_e - C_s
Cveg <- t (Cveg)
Cveg <- Cveg [, length (Cveg [1, ]):1]
Cveg <- as.matrix (Cveg)
Cveg [is.na (forest_mask_NA)] <- NA

# Read averaged MODIS and BETHYS NPP layers
#----------------------------------------------------------------------------------------
NPP <- read.table (file = 'NPP.txt', header = F)
NPP [NPP == -99] <- NA
NPP <- t (NPP)
NPP <- NPP [, length (NPP [1, ]):1]
NPP <- as.matrix (NPP)
NPP [is.na (forest_mask_NA)] <- NA

# Calculate tau 
#----------------------------------------------------------------------------------------
nta <- Cveg / NPP
nta_adj1 <- Cveg / (NPP * (forest_mask_NA / 100.0))
nta_adj2 <- Cveg / (NPP * (canopy_mask_NA / 100.0))

# Clean up
#----------------------------------------------------------------------------------------
rm (nc_name, ncin, threshold)
#========================================================================================