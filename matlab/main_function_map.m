% Example script to run CATMIP algorithm to sample a 2*Ncomponent-by-1 vector
% whose first Ncomponents elements are sizes, and second Ncomponent
% elements are corresponding abundances.
% Prior is uniform distribution on size, Dirichlet distribution on abundance.
% Data likelihood is multivariate normal (MVN) distribution
%
% Input parameters are highlighted with "***" in comments
%
% Please cite:
% Lapotre, M.G.A. et al. (2017), A probabilistic approach to remote
% compositional analysis of planetary surfaces, JGR: Planets, doi:
% 10.1002/2016JE005248
%
% Edits by (tc): Added support for determining the optimal subset of
% endmembers among a larger set, and saved result variables for analysis.
function [MAP_grain_sizes_in_microns, MAP_abundances_in_percent, MAP_endmembers_indices, R2, rmse, lt, wt, wmi] = ...
    main_function_map(results_dir, run_number, run_col, min_endmembers, max_endmembers, data, mixed_spectrum)

% Create the results directories for all of the outputs.
fit_dir = strcat(results_dir, 'BestFit_and_Residual\');
cov_dir = strcat(results_dir, 'Distributions_and_Correlations\');
spec_dir = strcat(results_dir, 'Spectra\');
var_dir = strcat(results_dir, 'Variables\');
mkdir(fit_dir);
mkdir(cov_dir);
mkdir(spec_dir);
mkdir(var_dir);

% Read in the endmember optical constants and all other parameters needed
% for the Hapke-based spectral unmixing.
D = load_data_map(data, mixed_spectrum);
min_endmembers = str2double(min_endmembers); % Minimum number of endmembers.
D.Ncomponents = str2double(max_endmembers); % Maximum number of endmembers.
[~, num_lib] = size(D.mineral_classes);

% Assemble the matrix of randomly-determined samples for the prior PDF.
% Each column represents one sample.
THETA = zeros(3*D.Ncomponents, D.Nm);
num_endmembers = randi([min_endmembers, D.Ncomponents], 1, D.Nm);

% First, pull the abundances from a Dirichlet distribution.
% Next, (pseudo)-randomly choose the endmembers.
% Finally, choose the grain sizes from a continuous uniform distribution,
% where the minimum and maximum values are endmember-specific.
for i = 1:D.Nm
    THETA(D.Ncomponents+1:D.Ncomponents+num_endmembers(i), i) = dirichletrnd_exc(D.prior_alpha(1:num_endmembers(i)), 1);
    THETA(2*D.Ncomponents+1:2*D.Ncomponents+num_endmembers(i), i) = sort(randperm(num_lib, num_endmembers(i)))';
    for j = 1:num_endmembers(i)
        endmember_index = THETA(2*D.Ncomponents+j, i);
        THETA(j, i) = unifrnd(D.prior_size_low(endmember_index), D.prior_size_high(endmember_index), 1, 1);
    end
end

% Function handle to evaluate prior PDF: accept model vector y and structure
% return ln p(theta)
llk_prior_handle = @llk_prior;

% Function handle to evaluate data likelihood: accept model vector y,
% structure, and ln p(theta), return ln p(D|theta) and (optionally) ln p(theta)
for_model = @likelihood;

% Setup complete

%%%%%%%%%%%% Inversion %%%%%%%%%%%%%%%
% Run CATMIP

tic
[THETA, LLK, nAccept, nReject] = ...
    catmip_dirichlet(D.Nm, D.Nsteps, THETA, D, llk_prior_handle, for_model);
toc

% Return values:
% THETA: Nparam x N x M+1 array. M is total number of cooling steps.
% THETA(:,:,1) are samples of the prior. THETA(:,:,end) is final
% solution. Each column of matrix at each cooling step is a single model.
% LLK: 3 x N x M+1 array. Same layout as THETA except each column of
% matrix at each cooling step are the three log-likelihoods: 
% ln [p(theta|D) p(D|theta) p(theta)]'
% nAccept: M+1 x 1 array containing number of accepted sampls
% nReject: M+1 x 1 array containing number of rejected sampls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save and print the MAP model (grain sizes, abundances, and endmembers).
lend = LLK(2, :, end);
ind = find(lend == max(lend));
if length(ind) > 1
    ind = ind(1);
end
MODS = THETA(:, :, end);
bestmod = MODS(:, ind);
MAP_grain_sizes_in_microns = bestmod(1:length(bestmod)/3)
MAP_abundances_in_percent = bestmod(length(bestmod)/3+1:length(bestmod)/3*2)
MAP_endmembers_indices = bestmod(length(bestmod)/3*2+1:end);
MAP_endmembers_labels = strings(length(bestmod)/3, 1);

MAP_grain_sizes_in_microns = MAP_grain_sizes_in_microns(MAP_grain_sizes_in_microns ~= 0)
MAP_abundances_in_percent = MAP_abundances_in_percent(MAP_abundances_in_percent ~= 0)
MAP_endmembers_indices = MAP_endmembers_indices(MAP_endmembers_indices ~= 0)
MAP_endmembers_labels = strings(length(MAP_endmembers_indices), 1);

for i = 1:length(MAP_endmembers_indices)
    index = MAP_endmembers_indices(i);
    if index > 0
        MAP_endmembers_labels(i) = D.labels(1, index);
    end
end
MAP_endmembers = strcat(string(MAP_endmembers_indices), ' (', MAP_endmembers_labels, ')')

close all; % Prevent too many figures from accumulating.

% Plot the MAP model, as well as all other models.
unique_mods = unique(MODS.', 'rows').';
[unique_mods_rows, ~] = size(unique_mods);
all_ab = unique_mods(unique_mods_rows/3+1:unique_mods_rows/3*2,:);
all_d = unique_mods(1:unique_mods_rows/3,:);
all_ind = unique_mods(unique_mods_rows/3*2+1:end,:);
[R2, rmse, lt, wt, wmi, wifi, f] = ...
        genspectralplot(MAP_abundances_in_percent, MAP_grain_sizes_in_microns, MAP_endmembers_indices, all_ab, all_d, all_ind, D.R_vs_SSA, D);

% Create images of each spectral plot, for convenient viewing.
if R2 >= 0
    saveas(gcf, strcat(fit_dir, num2str(R2, '%.6f'), '_', num2str(run_number), '_', num2str(run_col), '_BestFit_and_Residual.jpg'))
else
    saveas(gcf, strcat(fit_dir, '_', num2str(run_number), '_', num2str(run_col), '_BestFit_and_Residual.jpg'))
end

% Output the data and model spectra in ASCII format.
if R2 >= 0
    fileID = fopen(strcat(spec_dir, num2str(R2, '%.6f'), '_', num2str(run_number), '_', num2str(run_col), '_Spectra.txt'), 'wt');
else
    fileID = fopen(strcat(spec_dir, '_', num2str(run_number), '_', num2str(run_col), '_Spectra.txt'), 'wt');
end
[num_rows, ~] = size(lt);
num_endmembers = length(MAP_endmembers_indices);
fprintf(fileID, 'ENVI ASCII Plot File\n');
fprintf(fileID, 'Column 1: Wavelength\n');
fprintf(fileID, 'Column 2: Data\n');
fprintf(fileID, 'Column 3: Model\n');
for i = 1:num_endmembers
    index = MAP_endmembers_indices(i);
    fprintf(fileID, 'Column %d: %f %s\n', (i + 3), MAP_grain_sizes_in_microns(i), erase(D.labels(index), '\'));
end
for i = 1:num_endmembers
    index = MAP_endmembers_indices(i);
    fprintf(fileID, 'Column %d: %f %f %f %s\n', (i + 3 + num_endmembers), f(i), MAP_abundances_in_percent(i), MAP_grain_sizes_in_microns(i), erase(D.labels(index), '\'));
end
for i = 1:num_rows
    fprintf(fileID, '%f\t%f\t%f', lt(i)*1e6, wt(i), wmi(i)); % Write the data and model SSAs.
    for j = 1:num_endmembers
        fprintf(fileID, '\t%f', wifi(i, j)); % Write the component SSAs.
    end
    for j = 1:num_endmembers
        fprintf(fileID, '\t%f', f(j) * wifi(i, j)); % Write the component SSAs weighted by fractional cross section.
    end
    fprintf(fileID, '\n');
end
fclose(fileID);

%%%%% Todo: Don't know how to handle this yet for varying numbers of endmembers.
% gencorrelplot(MODS, bestmod, D);
% if R2 >= 0
%     saveas(gcf, strcat(cov_dir, num2str(R2, '%.6f'), '_', num2str(run_number), '_', num2str(run_col), '_Distributions_and_Correlations.jpg'))
% else
%     saveas(gcf, strcat(cov_dir, '_', num2str(run_number), '_', num2str(run_col), '_Distributions_and_Correlations.jpg'))
% end

% Save all of the relevant variables to reconstruct the results, along with confidence intervals.
if R2 >= 0
    save(strcat(var_dir, num2str(R2, '%.6f'), '_', num2str(run_number), '_', num2str(run_col), '_Results.mat'), 'MODS', 'bestmod', 'D')
else
    save(strcat(var_dir, '_', num2str(run_number), '_', num2str(run_col), '_Results.mat'), 'MODS', 'bestmod', 'D')
end

end
