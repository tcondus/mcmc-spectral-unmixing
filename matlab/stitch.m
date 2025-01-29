% Combine piecewise abundance and grain size maps into a single map.
function stitch

% Select two or more files to stitch together, and initialize the output.
[filename, pathname, ~] = uigetfile('*', 'Select Two or More Files', 'MultiSelect', 'on');
file = fullfile(pathname, filename{1, 1});
[data, hdr] = read_envi_data(file);
stitched_cube = data;

num_files = length(filename);
for i = 2:num_files
    file = fullfile(pathname, filename{1, i});
    [data, hdr] = read_envi_data(file);
    [data_rows, data_cols, data_layers] = size(data);
    for j = 1:data_rows
        for k = 1:data_cols
            pixel = data(j, k, :);
            if sum(pixel) > 0
                stitched_cube(j, k, :) = pixel;
            end
        end
    end
end

fprintf('Writing stitched cube file...');
current_time = now;
current_time_str = strcat(datestr(current_time, 'yyyy'), datestr(current_time, 'mm'), datestr(current_time, 'dd'), '_', ...
    datestr(current_time, 'HH'), datestr(current_time, 'MM'), datestr(current_time, 'SS'), datestr(current_time, 'FFF'));
results_filename = strcat(pathname, 'stitched_cube_', current_time_str);
write_envi_data(stitched_cube, results_filename, 'bsq');
fileID = fopen(strcat(results_filename, '.hdr'), 'at');
fprintf(fileID, 'map info = {%s}\n', hdr.map_info);
fprintf(fileID, 'projection info = {%s}\n', hdr.projection_info);
fprintf(fileID, 'coordinate system string = {%s}\n', hdr.coordinate_system_string);
fprintf(fileID, 'data ignore value = 0.00000000000000e+00\n');
if isfield(hdr, 'band_names')
    fprintf(fileID, 'band names = {\n');
    [~, num_bands] = size(hdr.band_names);
    for i = 1:num_bands-1
        fprintf(fileID, ' %s, \n', hdr.band_names{i});
    end
    fprintf(fileID, ' %s\n', hdr.band_names{i+1});
    fprintf(fileID, '}\n');
end
if isfield(hdr, 'wavelength')
    fprintf(fileID, 'wavelength = {\n');
    [~, num_bands] = size(hdr.wavelength);
    for i = 1:num_bands-1
        fprintf(fileID, ' %f, \n', hdr.wavelength(i));
    end
    fprintf(fileID, ' %f\n', hdr.wavelength(i+1));
    fprintf(fileID, '}\n');
end
fclose(fileID);
fprintf(' done.\n');

end
