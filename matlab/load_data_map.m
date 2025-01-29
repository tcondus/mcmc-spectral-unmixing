% Used by mcmc_map.m: Nm = 2000 and Nsteps = 10.
% Todo: Allow the user to specify these values in the run file instead.
function D = load_data_map(data, mixed_spectrum)

[data_rows, ~] = size(data);
cells = cell(data_rows, 5);
for i = 1:data_rows
    split = strsplit(data(i));
    cells{i, 1} = str2double(split(1));
    cells{i, 2} = str2double(split(2));
    cells{i, 3} = str2double(split(3));
    cells{i, 4} = split(4);
    cells{i, 5} = str2double(split(5));
end
[cell_rows, ~] = size(cells);

%-------------------------------------------------------------------------
% INPUT PARAMETERS
%-------------------------------------------------------------------------

% MINERAL ENDMEMBERS

mineral_classes = zeros(1, cell_rows);
for i = 1:cell_rows
    mineral_classes(i) = cells{i, 5};
end

% *** Populate the density vector.
rho = zeros(cell_rows, 1);
for i = 1:cell_rows
    rho(i) = cells{i, 1};
end
% *** Load optical constants for each laboratory endmember; these are three
% columns matrices (wavelengths in meters, n, and k)
for i = 1:cell_rows
    P{i} = dlmread(cells{i, 4});
    [~, labels(i), ~] = fileparts(cells{i, 4});
end
N = length(P);
% Sort them into n (real) and k (imaginary) refractive indices
for i = 1:N
    O = P{i};
    l = O(:, 1);
    nn = O(:, 2);
    kk = O(:, 3);
    lambda{i} = l;
    n{i} = nn;
    k{i} = kk;
    clear O l nn kk
end

% MCMC PARAMETERS
Nm = 2000;      % *** Length of each Markov chain;
Nsteps = 10;   % *** Number of samples per cooling step
prior_size_low = zeros(1, cell_rows);
prior_size_high = zeros(1, cell_rows);
for i = 1:cell_rows
    prior_size_low(1, i) = cells{i, 2};
    prior_size_high(1, i) = cells{i, 3};
end
prior_alpha = ones(1, length(rho)); % Alpha parameter for each component; taken as uniform as descibed in Lapotre et al. (2017); can be modified if necessary
% Note that currently proposal PDF draws from Dirichlet distribution
% using same alpha_i as the prior PDF.  This can be edited if necessary.
% (Dirichlet distribution: http://en.wikipedia.org/wiki/Dirichlet_distribution)
[~, Ncomponents] = size(unique(mineral_classes));

% SPECTRAL DATA
R1_vs_SSA0 = 0; % *** If reflectance data, set to 1; if SSA, set to 0;
if R1_vs_SSA0
    inc = 30; % *** incident angle for data reflectance spectrum
    e = 0; % *** emission angle for data reflectance spectrum
else
    inc = NaN;
    e = NaN;
end
% *** Data (two columns matrix with wavelength and reflectance/SSA)
S = dlmread(mixed_spectrum);
S(:, 1) = S(:, 1) .* 1e-6; % Convert wavelengths from microns to meters.
lam_SPEC = S(:, 1);
R_SPEC = S(:, 2);
if min(lam_SPEC) > 1
    lam_SPEC = lam_SPEC .* 1e-9;
end

% FIT PARAMETERS
Cov = 0.0005; % *** Covariance matrix = (Cd + Cp); increase to allow for more variance
min_lam = 0.75 .* 1e-6; % Minimum wavelength of range over which spectral fit is performed, in meters
max_lam = 2.5 .* 1e-6; % Maximum wavelength of range over which spectral fit is performed, in meters
% Note: Data and endmember spectra needed over this wavelength range

% Creating Data Structure
D.inc = inc;
D.e = e;
D.rho = rho';
D.lambda = lambda;
D.n = n;
D.k = k;
D.lam_SPEC = lam_SPEC;
D.R_SPEC = R_SPEC;
D.rho = rho'; 
D.Cinv = 1/Cov; 
D.Ncomponents = Ncomponents; 
D.prior_size_low = prior_size_low; 
D.prior_size_high = prior_size_high; 
D.prior_alpha = prior_alpha; 
D.Nm = Nm;
D.Nsteps = Nsteps;
D.R_vs_SSA = R1_vs_SSA0;
D.lmin = min_lam;
D.lmax = max_lam;
D.labels = strrep(labels, '_', '\_'); % Print the underscore in the legend.
D.mineral_classes = mineral_classes;

end
