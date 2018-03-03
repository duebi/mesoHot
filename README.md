# mesoHot

This repository contains MATLAB routines to analyze how changes in sea surface height are associated with changes in the biogeochemistry of Station ALOHA as measured during the Hawaii Ocean Time-series (HOT) for the 1993-2015 period.
Sea surface height was produced and distributed by the Copernicus Marine and Environment Monitoring Service (CMEMS, http://www.marine.copernicus.eu). HOT observations were downloaded from the Hawaii Ocean Time-series Data Organization & Graphical System (HOT-DOGS, http://hahana.soest.hawaii.edu/hot/hot-dogs/).

Requirements:
- Matlab ver 2013b or above
- Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 available rom http://www.TEOS-10.org
- lsqfitgm.m, lsqfitx.m and lsqfity.m available from http://www.mbari.org/index-of-downloadable-files/ (by E.T. Peltzer) 
- fdr_bh.m for the Benjamini & Hochberg procedure available from https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr-bh (v2.3 by David M. Groppe)
- readHotDogs.m, plot_REG_phase_h.m, plot_REG_results_h.m, plot_REG_table_h.m, redbluecmap.m, subplot_labels.m, zerocmap.m 

The auxiliary/ folder contains routines used within mesoHot.m and mesoHot_AR.m
The data/ folder contains altimetric data and HOT hydrographic and biogeochemical observations. This folder needs to be added to your MATLAB path before calling mesoHot.m or other routines.

Example of use (in MATLAB’s Command Window):

[rs,pvals,rs_AR,pvals_AR,pvals_rks] =  mesoHot('pro’,’Spearman’)

Produces a plot showing the extreme quartile analysis for Prochlorochoccus in the upper 175 m of the water column & the correlation coefficients from the linear analysis and the autoregression residuals analysis. Function outputs are the correllation coefficients and p values.

The script tables_mesoHot.m contains (in the first cell) the routine to produce the synthetic results from the correlations of multiple biogeochemical variables and multiple depths. The script mesoHOT_isopycnal.m is divided in different cells that do several isopycnal analyses.

Benedetto Barone and Ashley Coenen - Feb 2018
