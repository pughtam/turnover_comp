Scripts underlying the drought mortality evaluation in Table S5 of:
Thomas A. M. Pugh, Tim Rademacher, Sarah L. Shafer, Jörg Steinkamp, Jonathan Barichivich, Brian Beckage, Vanessa Haverd, Anna Harper, Jens Heinke, Kazuya Nishina, Anja Rammig, Hisashi Sato, Almut Arneth, Stijn Hantson, Thomas Hickler, Markus Kautz, Benjamin Quesada, Benjamin Smith, Kirsten Thonicke. Understanding the uncertainty in global forest carbon turnover.

For queries contact Jörg Steinkamp, joerg.steinkamp@uni-mainz.de (author) or Thomas A. M. Pugh, t.a.m.pugh@bham.ac.uk

These scripts use the EcoTools package which can be downloaded from https://gitlab.com/jsteinkamp/EcoTools. Data on drought mortality event locations is located in "Data", Scripts are located in "R".

Instructions:
1.) run 'R/read_locations.R', which reads the Excel table and saves it as 'data/dm_events.RData' ("drought mortality events")
2.) run 'R/globalmap.R' to see location of drought events (not strictly needed)
3.) run 'R/extract_mort.R', which extracts the timeseries from each individual model and averages it over the whole region.
4.) run R/visualize.R'
