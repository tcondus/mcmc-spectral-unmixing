% Wrapper code to run MCMC for every pixel in a CRISM scene, creating
% abundance and grain size maps. The set of endmembers used will be
% determined by the algorithm.
function mcmc_map

% The user may adjust these parameters to tweak performance.
R2_threshold = 0.9900; % The minimum R2 score to accept.
max_trials = 1; % How many runs to perform until the minimum R2 is met.
%max_trials = 3;

% Select a run file containing (in row order):
% (1) The prefix that will identify this run.
% (2) The absolute path to the results directory.
% (3) The absolute path to the CRISM scene.
% (4) The minimum number of endmembers to pull from the spectral library.
% (5) The maximum number of endmembers to pull from the spectral library.
% (6) The minimum wavelength of the spectral range.
% (7) The maximum wavelength of the spectral range.
% (8) The minimum row number of the data cube.
% (9) The maximum row number of the data cube.
% (10) The minimum column number of the data cube.
% (11) The maximum column number of the data cube.
% (12-end) The absolute path(s) to the endmember class files. For now, only
%          a single file containing the endmember information is supported.
[filename, pathname, ~] = uigetfile({'*.txt', 'Text Files'; '*.*', 'All Files'}, 'Load run file');
if filename == 0
    disp('Missing run file.');
    return;
end
run_fileID = fopen(fullfile(pathname, filename), 'r');
run_file_data = textscan(run_fileID, '%s');
[num_lines, ~] = size(run_file_data{1, 1});
fclose(run_fileID);

% Hard-coded run file structure.
run_prefix = run_file_data{1, 1}{1, 1};
results_path = run_file_data{1, 1}{2, 1};
ssa_join_warp = run_file_data{1, 1}{3, 1};
min_endmembers = run_file_data{1, 1}{4, 1};
max_endmembers = run_file_data{1, 1}{5, 1};
data_min = str2double(run_file_data{1, 1}{6, 1});
data_max = str2double(run_file_data{1, 1}{7, 1});
row_min = str2double(run_file_data{1, 1}{8, 1});
row_max = str2double(run_file_data{1, 1}{9, 1});
col_min = str2double(run_file_data{1, 1}{10, 1});
col_max = str2double(run_file_data{1, 1}{11, 1});
grid_heading = [];
for i = 12:num_lines % Iterate over any number of endmember classes.
    endmember_class = run_file_data{1, 1}{i, 1};
    endmember_class_fileID = fopen(endmember_class, 'r');
    endmember_data = textscan(endmember_class_fileID, '%s %s %s %s');
    fclose(endmember_class_fileID);
    [num_lib, ~] = size(endmember_data{1, 1});
    endmember_array = strings(num_lib, 1);
    for j = 1:num_lib
        endmember_array(j) = strcat(num2str(endmember_data{1, 1}{j, 1}), " ", ...
                num2str(endmember_data{1, 2}{j, 1}), " ", ...
                num2str(endmember_data{1, 3}{j, 1}), " ", ...
                endmember_data{1, 4}{j, 1}, " ", num2str(i - 5));
    end
    grid_heading = [grid_heading, endmember_array'];
end
[~, total_endmembers] = size(grid_heading);

% Prepare output directory.
current_time = now;
current_time_str = strcat(datestr(current_time, 'yyyy'), datestr(current_time, 'mm'), datestr(current_time, 'dd'), '_', ...
    datestr(current_time, 'HH'), datestr(current_time, 'MM'), datestr(current_time, 'SS'), datestr(current_time, 'FFF'));
run_prefix_unique = strcat(run_prefix, '_', current_time_str);
results_dir = strcat(results_path, run_prefix_unique, '\');
mkdir(results_dir);

% Read in the CRISM data.
[data, hdr] = read_envi_data(ssa_join_warp);
[data_rows, data_cols, ~] = size(data);

% Number of backplanes: first n for the abundances, next n for the grain sizes, and final 2 for R2 and RMSE.
num_backplanes = (total_endmembers * 2) + 2;
results_cube = zeros(data_rows, data_cols, num_backplanes);
num_model_bands = 0;

% Run the MCMC algorithm for each pixel in the CRISM scene subset.
for i = row_min:row_max
    for j = col_min:col_max
        temp_filename = strcat(results_dir, run_prefix_unique, '_', num2str(i), '_', num2str(j), '.txt');
        tempID = fopen(temp_filename, 'wt');
        ssa = data(i, j, :);
        wav = hdr.wavelength;
        [~, num_wav] = size(wav);
        for k = 1:num_wav
            current_wav = wav(1, k);
            current_ssa = ssa(1, 1, k);
            if (current_wav > data_min && current_wav < data_max)
                fprintf(tempID, '%f\t%f\n', current_wav, current_ssa);
            end
        end
        fclose(tempID);
        
        % Only consider the pixels that have valid data.
        if sum(ssa) ~= 0
            current_run = grid_heading';
            current_run = current_run(~cellfun('isempty', current_run));
            try
                highest_R2 = 0;
                highest_R2_index = 0;
                trials = 0;
                while true
                    % Run the main MCMC function.
                    trials = trials + 1;
                    [grain_sizes, abundances, endmembers, R2, rmse, lt, wt, wmi] = ...
                            main_function_map(results_dir, i, j, min_endmembers, max_endmembers, current_run, temp_filename);
                    if num_model_bands == 0
                        num_model_bands = length(wmi);
                        model_cube = zeros(data_rows, data_cols, num_model_bands);
                        model_cube_trials = zeros(data_rows, data_cols, num_model_bands, 5);
                        results_cube_trials = zeros(data_rows, data_cols, num_backplanes, 5);
                    end

                    % Keep track of the run with the highest R2 score among
                    % all trials.
                    if highest_R2 < R2
                        highest_R2 = R2;
                        highest_R2_index = trials;
                    end
                    %if lowest_rmse <= 0.0020
                    if highest_R2 >= R2_threshold
                        for k = 1:length(endmembers)
                            results_cube(i, j, endmembers(k)) = abundances(k);
                            results_cube(i, j, total_endmembers + endmembers(k)) = grain_sizes(k);
                        end
                        results_cube(i, j, num_backplanes-1) = R2;
                        results_cube(i, j, num_backplanes) = rmse;
                        model_cube(i, j, :) = wmi;
                        break;
                    else
                        for k = 1:length(endmembers)
                            results_cube_trials(i, j, endmembers(k), trials) = abundances(k);
                            results_cube_trials(i, j, total_endmembers + endmembers(k), trials) = grain_sizes(k);
                        end
                        results_cube_trials(i, j, num_backplanes-1, trials) = R2;
                        results_cube_trials(i, j, num_backplanes, trials) = rmse;
                        model_cube_trials(i, j, :, trials) = wmi;
                    end
                    
                    % If none of the trials meets the R2 score threshold,
                    % then choose the run with the highest score.
                    if trials == max_trials
                        results_cube(i, j, :) = results_cube_trials(i, j, :, highest_R2_index);
                        model_cube(i, j, :) = model_cube_trials(i, j, :, highest_R2_index);
                        break;
                    end
                end
                
                fprintf('Run (%d, %d) complete.\n', i, j);
            catch exception
                fprintf('Error: Run (%d, %d) incomplete.\n', i, j);
                disp(getReport(exception, 'extended'));
            end
        else
            fprintf('Run (%d, %d) complete, no data.\n', i, j);
        end
    end
end

% Write the abundance and grain size results in ENVI format.
fprintf('Writing abundance and grain size results file...');
results_filename = strcat(results_dir, run_prefix_unique, '_results');
write_envi_data(results_cube, results_filename, 'bsq');
fileID = fopen(strcat(results_filename, '.hdr'), 'at');
fprintf(fileID, 'map info = {%s}\n', hdr.map_info);
fprintf(fileID, 'projection info = {%s}\n', hdr.projection_info);
fprintf(fileID, 'coordinate system string = {%s}\n', hdr.coordinate_system_string);
fprintf(fileID, 'data ignore value = 0.00000000000000e+00\n');
fprintf(fileID, 'band names = {\n');
for i = 1:total_endmembers
    [~, endmember, ~] = fileparts(grid_heading(i));
    fprintf(fileID, ' %s abundance, \n', endmember);
end
for i = 1:total_endmembers
    [~, endmember, ~] = fileparts(grid_heading(i));
    fprintf(fileID, ' %s grain size, \n', endmember);
end
fprintf(fileID, ' R^2, \n');
fprintf(fileID, ' RMSE\n');
fprintf(fileID, '}\n');
fclose(fileID);
fprintf(' done.\n');

% Write the model spectra cube in ENVI format.
fprintf('Writing model spectra cube file...');
model_filename = strcat(results_dir, run_prefix_unique, '_model');
write_envi_data(model_cube, model_filename, 'bsq');
fileID = fopen(strcat(model_filename, '.hdr'), 'at');
fprintf(fileID, 'map info = {%s}\n', hdr.map_info);
fprintf(fileID, 'projection info = {%s}\n', hdr.projection_info);
fprintf(fileID, 'coordinate system string = {%s}\n', hdr.coordinate_system_string);
fprintf(fileID, 'data ignore value = 0.00000000000000e+00\n');
fprintf(fileID, 'wavelength = {\n');
for i = 1:length(lt)-1
    fprintf(fileID, ' %f,\n', lt(i)*1e6);
end
fprintf(fileID, ' %f\n', lt(end)*1e6);
fprintf(fileID, '}\n');
fclose(fileID);
fprintf(' done.\n');

end
