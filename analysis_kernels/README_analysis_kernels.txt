Scripts underlying kernel figures (Fig. 2 and 5) and preprocessing of satellite-based estimates in:
Thomas A. M. Pugh, Tim Rademacher, Sarah L. Shafer, Jörg Steinkamp, Jonathan Barichivich, Brian Beckage, Vanessa Haverd, Anna Harper, Jens Heinke, Kazuya Nishina, Anja Rammig, Hisashi Sato, Almut Arneth, Stijn Hantson, Thomas Hickler, Markus Kautz, Benjamin Quesada, Benjamin Smith, Kirsten Thonicke. Understanding the uncertainty in global forest carbon turnover.

For queries, contact Tim T. Rademacher, trademacher@fas.harvard.edu or Thomas A. M. Pugh, t.a.m.pugh@bham.ac.uk

The scripts are divided into scripts that read and wrangle the data for each TBM (e.g. r_CABLE.R and CABLE_tau.R) and a script that plots the density kernels (plot_density_kernels_Fig3_and_Fig6.R). The plotting script sources the read and wrangle scripts and thus relies on them.

read_and_plot_observational_data.R processes the satellite-based estimates.

03.12.19
