Scripts underling the majority of the figures in the analysis.
Scripts were developed for Matlab 2018a by T. A. M. Pugh.

Basic data processing to read in model data from netcdfs, calculate global statistics, apply masks and save to *.mat binary files
- turnover_pool_flux_read.m

Helper functions for data processing. See file headers for details.
- global_grid_area.m
- esa_05_landcover.m
- esa_jules_landcover.m
- get_forest_type.m
- get_stocks_fluxes.m
- get_closed_can_mask.m

Additional dependencies required from elsewhere
- hansen_forest_frac_0p5deg.nc4 (available from doi: 10.18161/disturbance_forestmask.201905)

The following scripts are used to produce the indicated plots

global_turn_table.m
global_obs_totals.m
Table 3

global_turn_frac_map.m
Fig. 2
Fig. S3
Fig. S6

global_turn_frac_sd_space_time_comb.m
global_turn_frac_sd.m
global_turn_frac_sd_time.m
Fig. 4
Fig. S8

turn_breakdown_bar.m
Fig.5

dom_mort_mech_map.m
Fig. 7

ipsl_line_plots.m
Fig. 8

global_turn_diff_map_ipsl.m
Fig. 9

forest_cover_comparison.m
Fig. S1

aleaf_vs_SLA_tradeoffs.m
Fig. S4

lai_map_plot.m
Fig. S5

phen_map_plot.m
Fig. S7

tropical_plot_comparison.m
Fig. S9

phen_map_plot_ipsl.m
Fig. S10

biom_flux_turn_time.m
Fig. S11-S17

29.11.19
