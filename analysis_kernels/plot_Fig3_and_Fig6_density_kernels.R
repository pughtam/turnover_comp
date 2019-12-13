#========================================================================================
# Plot figures 3 and 6 for Pugh et al. (2019) titled 
# "Understanding the uncertainty in global forest carbon turnover".
#
# This code was writen by Tim T. Rademacher (trademacher@fas.harvard.edu).
#
# Last modified: 2019/12/06
#----------------------------------------------------------------------------------------

# Load observational data
#----------------------------------------------------------------------------------------
source ('read_and_plot_observational_data.R')

# Set colours for phenologies
#----------------------------------------------------------------------------------------
colours <- c ('#91b9a4', '#106470', '#EB99A9', '#55A51C', '#EA7125', '#8F2BBC', '#D6083B')

# phenologies
#----------------------------------------------------------------------------------------
# 1 = needleleaved evergreen
# 2 = needleleaved deciduous
# 3 = boreal broadleaved deciduous
# 4 = temperate broadleaved evergreen
# 5 = temperate broadleaved deciduous
# 6 = tropical broadleaved evergreen
# 7 = tropical broadleaved raingreen

# Adjust categories for observational phenology mask
#----------------------------------------------------------------------------------------
pheno_mask_obs <- matrix (as.numeric (rep (NA, 720 * 360)), nrow = 720)
pheno_mask_obs [pheno_mask == 6] <- 1 # NE
pheno_mask_obs [pheno_mask == 7] <- 2 # ND
pheno_mask_obs [pheno_mask == 8] <- 3 # boreal BD or set it to 3
pheno_mask_obs [pheno_mask == 4] <- 4 # temperate BE
pheno_mask_obs [pheno_mask == 5] <- 5 # temperate BD
pheno_mask_obs [pheno_mask == 1] <- 6 # tropical BE
pheno_mask_obs [pheno_mask == 2] <- 7 # tropical BD
pheno_mask_obs [pheno_mask == 9] <- NA # tropical BD
  
# Function to add polygon
#----------------------------------------------------------------------------------------
add.polygon <- function (x1, y1, y2, i = 1) {
  polygon (x   = c (x1, rev (x1)),
           y   = c (y1, rev (y2)),
           lty = 0,
           col = colours [i])
  
}

# Function to add the mean 
#----------------------------------------------------------------------------------------
add.mean <- function (i, nta, pheno_mask) {
  points (x = mean (nta [pheno_mask == i], na.rm = T),
          y = -0.01,
          lwd = 2, pch = 19,
          col = colours [i])
}

# Function to plot density kernel by phenology for a TBM
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
        xlab = expression (paste (tau[mort],' (years)')),
        ylab = 'density',
        col = colours [1],
        lwd = 0,
        main = '',
        xlim = if (eq) {c (0, 35)} else  {c (0, 70)},
        ylim = if (eq) {c (-0.01, 0.18)} else {c (-0.01, 0.30)})
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
        y = if (eq) {0.17} else {0.29},
        labels = string,
        pos = 4)
  legend (x = if (eq) {28} else {55},
          y = if (eq) {0.18} else {0.30},
          legend = c ('NE','ND','boreal BD','temperate BE','temperate BD',
                      'tropical BE','tropical BR'),
          col = colours,
          pch = 19,
          cex = 0.8,
          box.lty = 0)
}

# Get CABLE fields
#----------------------------------------------------------------------------------------
source ('r_CABLE-POP.R')
rm (forest_mask_CABLE, forest_mask_CABLE_NA)

# Plot CABLE
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = CABLE_CRUN_nta,
              pheno_mask = CABLE_pheno_mask,
              string = '(a) CABLE')

# Get JULES fields
#----------------------------------------------------------------------------------------
source ('r_JULES.R')
rm (forest_mask_JULES, forest_mask_JULES_NA)

# Plot JULES
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = JULES_CRUN_nta,
              pheno_mask = JULES_pheno_mask,
              string = '(b) JULES')

# Get LPJ-GUESS fields
#----------------------------------------------------------------------------------------
source ('r_LPJ-GUESS.R')
rm (forest_mask_LPJ_GUESS, forest_mask_LPJ_GUESS_NA)

# Plot LPJ-GUESS
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = LPJ_GUESS_CRUN_nta,
              pheno_mask = LPJ_GUESS_pheno_mask,
              string = '(c) LPJ-GUESS')

# Get LPJmL fields
#----------------------------------------------------------------------------------------
source ('r_LPJmL.R')
rm (forest_mask_LPJmL, forest_mask_LPJmL_NA)

# Plot LPJmL
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = LPJmL_CRUN_nta,
              pheno_mask = LPJmL_pheno_mask,
              string = '(d) LPJmL')

# Get ORCHIDEE fields
#----------------------------------------------------------------------------------------
source ('r_ORCHIDEE.R')
rm (forest_mask_ORCHIDEE, forest_mask_ORCHIDEE_NA)

# Plot ORCHIDEE
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = ORCHIDEE_CRUN_nta_adj,
              pheno_mask = ORCHIDEE_pheno_mask,
              string = '(e) ORCHIDEE')

# Get SEIB fields
#----------------------------------------------------------------------------------------
source ('r_SEIB.R')
rm (forest_mask_SEIB, forest_mask_SEIB_NA)

# Plot SEIB
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = SEIB_CRUN_nta,
              pheno_mask = SEIB_pheno_mask,
              string = '(f) SEIB')

# Plot observational fields
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = nta_adj1,
              pheno_mask = pheno_mask_obs,
              string = '(g) Satellite-based')

# Get CABLE fields
#----------------------------------------------------------------------------------------
source ('CABLE_tau.R')

# Plot CABLE
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = CABLE_CRUN_tau,
              pheno_mask = CABLE_pheno_mask,
              string = 'CABLE',
              eq = FALSE)

# Get JULES fields
#----------------------------------------------------------------------------------------
source ('JULES_tau.R')

# Plot JULES
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = JULES_CRUN_tau,
              pheno_mask = JULES_pheno_mask,
              string = 'JULES',
              eq = F)

# Get LPJ-GUESS fields
#----------------------------------------------------------------------------------------
source ('LPJ-GUESS_tau.R')
 
# Plot LPJ-GUESS
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = LPJ_GUESS_CRUN_tau,
              pheno_mask = LPJ_GUESS_pheno_mask,
              string = 'LPJ-GUESS',
              eq = F)

# Get LPJmL fields
#----------------------------------------------------------------------------------------
source ('LPJmL_tau.R')

# Plot LPJmL
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = LPJmL_CRUN_tau,
              pheno_mask = LPJmL_pheno_mask,
              string = 'LPJmL',
              eq = F)
 
# Get ORCHIDEE fields
#----------------------------------------------------------------------------------------
source ('ORCHIDEE_tau.R')
 
# Plot ORCHIDEE
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = ORCHIDEE_CRUN_tau,
              pheno_mask = ORCHIDEE_pheno_mask,
              string = 'ORCHIDEE',
              eq = F)

# Get SEIB fields
#----------------------------------------------------------------------------------------
source ('SEIB_tau.R')

# Plot SEIB
#----------------------------------------------------------------------------------------
pl.nta.pheno (nta = SEIB_CRUN_tau,
              pheno_mask = SEIB_pheno_mask,
              string = 'SEIB',
              eq = F)

#========================================================================================