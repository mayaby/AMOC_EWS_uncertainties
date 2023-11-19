This is the github repository for the Nature Communications paper "Uncertainties in critical slowing down indicators of observation-based fingerprints of the Atlantic Overturning Circulation".

If you have any questions please email maya.ben-yami@tum.de

--------------------------

## Plotting notebooks

Two jupyter notebooks contain the python code for plotting the figures from the paper. Main_Plots.ipynb has the code for the main text figures, and SI_Plots.ipynb for the SI figures. 

## Analysis scripts

The folder analysis_scripts contains python and bash scripts which introduce all the code used to conduct the analysis in the paper. Note - most of these scripts are not operational, but rather introduce functions/ show examples of the code used. To get the results in the paper these scripts need to be modified to fit individual file structures, which differ across datasets. All the information needed to reproduce the analysis is contained in these scripts. 

* EWS_functions.py - functions for calculation CSD indicators (EWS)
* EWS_wrappers.py - wrappers for the above functions, with or without detrending
* count_profiles.py - goes through the monthly salinity profile files in EN4 and counts the number of "good" profiles in a given lat-lon gridcell in the analysis grid. Only used for visualisation in figure S5
* example_parallel_code.py - an example of the mpi4py parallel codes that were used to run most of the CSD analysis in the paper
* example_shift_grid.py - an example of the grid pre-processing necesseary for correct cross hatching in plots like figure 6
* generate_surrogates.py -  contains the wrapper functions to: generate fourier surrogates, generate AR2 surrogates and to modify a surrogate array using the salinity weights method described in the paper
* get_trends.py - calculates the linear trend of input timeseries using xarray
* process_files.sh - all the various cdo scripts used to pre-process files (select areas, get annual average, etc.)



