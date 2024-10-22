This code provides an analysis of hospital contact data, with the purpose of measuring the extent to which hospital contact rates are dependent on density. 

Part 1 of the analysis (files beginning with "1") conducts the statistical inference to estimate the parameter "phi" which determines the extent of density dependence. 
This has a main analysis ("1.1_clusterScript_multiTarget.R"), an analysis disaggregated by the categories of hospital staff ("1.2_clusterScript_catHosp.R") and a sensitivity analysis in which some key assumptions are relaxed ("1.3_clusterScript_sensAnal.R"). These are designed to be run in batch format, with a single integer argument (analysis_index) being provided to each script to select the relevant subset of the data for analysis. However, the outputs of these analysis are provided in the folder "analysis". 

The actual estimation is conducted in the script "1.0_MCMC_multiTarget.R", which is called by each of these three analysis scripts. 

Part 2 of the analysis combines and loads these analyses, then produces AIC tables to compare the different possible models of density dependence, and finally producing the figures for the paper, all of which go into the folder "output". 

Part 3 produces an additional figure using a toy model for an illustrative display of the effect of the density dependence assumption. 

