# mcmc-spectral-unmixing

## Introduction

Fork of the [Markov chain Monte Carlo (MCMC) spectral unmixing MATLAB code](https://drive.google.com/file/d/17QSbC3j323dn57jjdatoTlBnSUZmhKH5/view) described by [Lapôtre et al. (2017)](https://doi.org/10.1002/2016JE005248).

Similar to other unmixing algorithms, MCMC spectral unmixing decomposes a mixed spectrum into its constituent endmembers, each weighted by their fractional abundance. The key advantage of using MCMC is that it returns probability density functions of abundances for each endmember (i.e., multiple sets of best-fit parameters), as opposed to just a single answer. Therefore, this method addresses the issue of nonuniqueness of model solutions for the spectral unmixing problem.

A limitation of the original script is that it will attempt to fit all user-provided endmembers into the model solutions, which requires manual curation of the endmember library to adjust for poor fits. Additionally, the procedure can only unmix one spectrum per run.

This fork thus provides the following extensions to the MCMC spectral unmixing code:

- Adds wrapper functions for performing multiple unmixing runs. Useful for assessing the consistency of unmixing a single spectrum (`mcmc.m`), or unmixing numerous pixels across a spectral image cube (`mcmc_map.m`).
- Dynamically determine the optimal set of endmembers at runtime (sparse unmixing). The user supplies a large library of endmembers, and subsets will automatically be selected in search of best-fit models. Beneficial when unmixing image cubes that feature great spectral diversity.

## Running the program

Executing `mcmc.m` with an input file will run the unmixing algorithm 100 times for a single mixed spectrum. The outputs include best-fit spectral and error plots, ASCII files showing the constituent endmembers, and MAT-files containing saved variables related to the unmixing procedure (in case of further analysis). All of these results will be organized into appropriately-named subfolders. In addition, summary files listing the abundances and grain sizes of all contributing endmembers will be outputted to the specified home folder.

Executing `mcmc_map.m` with an input file will run the unmixing algorithm for each pixel in a given spectral image cube. The outputs for each pixel include best-fit spectral and error plots, ASCII files showing the constituent endmembers, and MAT-files containing saved variables related to the unmixing procedure (in case of further analysis). All of these results will be organized into appropriately-named subfolders. Two image cubes will also be outputted to the home folder: (1) the model cube containing the maximum a posteriori (MAP) best-fit spectra for each pixel, and (2) the results cube containing abundance and grain size distribution maps for each endmember. To save processing time, it is recommended that the user subset the image cube and perform parallel MATLAB runs on each chunk, rather than running the code across the entire cube at once. Then the user may combine the chunked results using `stitch.m`. (Todo: Find a more streamlined way to improve parallelism for large runs.)

Note: See the comments for a description of how to setup the input files. (Todo: Add template files to illustrate sample usage.)

## Scientific applications

MCMC spectral unmixing has been used in a number of publications for unmixing visible and near-infrared (VNIR) spectra of laboratory mixtures of minerals, as well as Mars Reconnaissance Orbiter CRISM data. This code should also support the unmixing of other (VNIR) datasets, such as Mars Express OMEGA orbiter data, and Mars 2020 SuperCam rover data. The unmixing of mid-infrared datasets like Mars Global Surveyor TES and Mars Exploration Rover Mini-TES should also be theoretically possible using this method, provided that the wavelength parameters in `load_data.m` and `load_data_map.m` are adjusted accordingly.

Note: This unmixing algorithm relies on the mixed spectrum being in units of single-scattering albedo (SSA), and the input endmember library in terms of their per-wavelength optical constants. This is because the underlying reflectance model is given by the [Hapke function](https://doi.org/10.1029/JB086iB04p03039), which relates SSAs, optical constants, and reflectances for endmember mixtures. (Todo: Add utility scripts for converting endmember spectra from reflectance to optical constants.)

## Future work

- Restore functionality for generating correlation plots (temporarily omitted because of uncertainty of implementation for sparse unmixing)
- Make output plots less hardcoded for wavelength ranges
- Find a more streamlined way to improve parallelism for large runs.
- Add template files to illustrate sample usage
- Add utility scripts for converting endmember spectra from reflectance to optical constants
- Add more robust error checking
- Allow the user to more easily adjust MCMC hyperparameters (which are currently hardcoded)
- Develop a graphical interface
- Port to Python for better accessibility

## Credits

The original MCMC spectral unmixing code was developed by Mathieu Lapôtre et al., and was modified and redistributed here with permission. ENVI-related helper scripts `parseParameters`, `read_envi_data`, `read_envi_header`, `write_envi_data`, and `write_envi_header` were originally written by Daniel Politte and used here with permission.
