% Wrapper code to perform 100 MCMC runs where the endmembers (and number of
% endmembers) will be determined by the algorithm. The number of endmembers
% in the final solution will fall within the min-max range specified by the
% user.
function mcmc

% Select a run file containing (in row order):
% (1) The prefix that will identify this run.
% (2) The absolute path to the results directory.
% (3) The absolute path to the data (or mixed) spectrum.
% (4) The minimum number of endmembers to pull from the spectral library.
% (5) The maximum number of endmembers to pull from the spectral library.
% (6-end) The absolute path(s) to the endmember class files. For now, only
%         a single file containing the endmember information is supported.
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
mixed_spectrum = run_file_data{1, 1}{3, 1};
min_endmembers = run_file_data{1, 1}{4, 1};
max_endmembers = run_file_data{1, 1}{5, 1};
%data_min = str2double(run_file_data{1, 1}{6, 1}); % Todo: Minimum wavelength for mixed spectrum.
%data_max = str2double(run_file_data{1, 1}{7, 1}); % Todo: Maximum wavelength for mixed spectrum.
grid_heading = [];
for i = 6:num_lines % Iterate over any number of endmember classes.
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

% Prepare output data structures.
n = 100; % Todo: Allow the user to specify the number of MCMC runs.
[~, total_endmembers] = size(grid_heading);
MAP_grain_sizes_grid = zeros(n, total_endmembers);
MAP_abundances_grid = zeros(n, total_endmembers);

% Prepare output files and directories.
current_time = now;
current_time_str = strcat(datestr(current_time, 'yyyy'), datestr(current_time, 'mm'), datestr(current_time, 'dd'), '_', ...
    datestr(current_time, 'HH'), datestr(current_time, 'MM'), datestr(current_time, 'SS'), datestr(current_time, 'FFF'));
run_prefix_unique = strcat(run_prefix, '_', current_time_str);
results_dir = strcat(results_path, run_prefix_unique, '\');
mkdir(results_dir);
fileID1 = fopen(strcat(results_dir, run_prefix_unique, '_abundances.txt'), 'wt');
fileID2 = fopen(strcat(results_dir, run_prefix_unique, '_grain_sizes.txt'), 'wt');
fprintf(fileID1, '%s\n', grid_heading(1:end));
fprintf(fileID2, '%s\n', grid_heading(1:end));

% Execute all runs.
for i = 1:n
    fprintf('Run %d starting...\n', i);
    current_run = grid_heading';
    current_run = current_run(~cellfun('isempty', current_run));
    fprintf(fileID1, '%d\t', i);
    fprintf(fileID2, '%d\t', i);
    try
        [grain_sizes, abundances, endmembers, R2, rmse] = ...
                main_function(results_dir, i, min_endmembers, max_endmembers, current_run, mixed_spectrum);
        N = length(endmembers);
        for j = 1:N
            index = endmembers(j);
            MAP_grain_sizes_grid(i, index) = grain_sizes(j);
            MAP_abundances_grid(i, index) = abundances(j);
        end
        fprintf(fileID1, '%f\t%f\t', R2, rmse);
        fprintf(fileID2, '%f\t%f\t', R2, rmse);
        for j = 1:total_endmembers
            fprintf(fileID1, '%f\t', MAP_abundances_grid(i, j));
            fprintf(fileID2, '%f\t', MAP_grain_sizes_grid(i, j));
        end
        fprintf('Run %d complete.\n', i);
    catch exception
        fprintf('Error: Run %d incomplete.\n', i);
        disp(getReport(exception, 'extended'));
    end
    fprintf(fileID1, '\n');
    fprintf(fileID2, '\n');
end
fclose(fileID1);
fclose(fileID2);

end
