function write_envi_header(filename,header_data)
% WRITE_ENVI_HEADER Generate an ENVI header with given filename & contents
% Creates an ENVI header file, _filename_, which contains only and exactly
% the data described in _header_data_ (with the single exception that a
% generic description is added if none is given in that structure).
% This function assumes that the header_data is consistant: e.g., that
% header_data.wavelength contains wavelengths in whichever unit is named in
% header_data.wavelength_units.

% There are some properties with 'Matlab' versions. We need to convert
% these values to ENVI-friendly values, and disregard the 'Matlab' versions
% when we write to the header.

if isfield(header_data, 'byte_order_matlab')
    switch header_data.byte_order_matlab
        case 'ieee-le'
            header_data.byte_order = 0;
        case 'ieee-be'
            header_data.byte_order = 1;
        otherwise
            error('The given ENVI byte order code, ''%s'', is invalid. ieee-le and ieee-be are the only valid values.', header_data.byte_order_matlab);
    end
    
    % Delete the Matlab-specific field so we won't accidentally write it
    % into the file
    header_data = rmfield(header_data, 'byte_order_matlab');
    
end

% If there's a Matlab data type listed, modify the ENVI data type number to
% match it. (If there isn't, we can just use the data_type value that's
% presumably there.)
if isfield(header_data, 'data_type_matlab')
    switch header_data.data_type_matlab
        case 'uint8'
            header_data.data_type = 1;
        case 'int16'
            header_data.data_type = 2;
        case 'int32'
            header_data.data_type = 3;
        case 'float'
            header_data.data_type = 4;
        case 'single'
            header_data.data_type = 4;
        case 'double'
            header_data.data_type = 5;
        case 'uint16'
            header_data.data_type = 12;
        case 'uint32'
            header_data.data_type = 13;
        case 'int64'
            header_data.data_type = 14;
        case 'uint64'
            header_data.data_type = 15;
        otherwise
            error('There is no know ENVI datatype code for the type ''%s''. If this is a valid type, please update this function.', header_data.data_type_matlab);
    end
    
    % Delete the Matlab-specific field so we won't accidentally write it
    % into the file
    header_data = rmfield(header_data, 'data_type_matlab');
    
end

% If there's no description, add a generic one

if ~isfield(header_data, 'description')
    header_data.description = 'ENVI header saved from Matlab';
end

% If there's a wavelength list, and the units are Nanometers, transform the
% units to be Micrometers instead (project policy)
if isfield(header_data, 'wavelength')
    if strcmpi(header_data.wavelength_units,'Nanometers')
        header_data.wavelength = header_data.wavelength / 1000; % nm -> um
        header_data.wavelength_units = 'Micrometers';
    end
end

% Any fields in our file that are mentioned in this list should be in this
% order, with any other fields coming after (and we don't care about their
% internal order)
% The ordering of fields we'll be using:
desired_order = {'description','samples','lines','bands','data_type','interleave','file_type', ...
    'sensor_type','header_offset','byte_order','default_bands','map_info','projection_info', ...
    'coordinate_system_string','band_names','wavelength_units','wavelength','fwhm','bbl', ...
    'data_ignore_value'};

header_fields = fieldnames(header_data);

header_fields = permute_strings_array(header_fields, desired_order);

fid = fopen(filename,'w');
fprintf(fid,'ENVI\r\n');

% Per field, write into the file
for i=1:length(header_fields)
    field_name = header_fields{i};
    field_name_for_writing = unclean_field_name(field_name);
    
    field_value = header_data.(field_name);
    
    % Determine what sort of value this is, and print it accordingly
    
    % preprocessing: if it's an array of numbers or boolean values, convert
    % to a cell array of strings
    if (isnumeric(field_value) || islogical(field_value)) && ~isscalar(field_value)
        field_value_strings = cell(1,length(field_value));
        for j=1:length(field_value)
            field_value_strings{j} = num2str(field_value(j));
        end
        field_value = field_value_strings;
    end
    
    
    if ischar(field_value) || (isnumeric(field_value) && isscalar(field_value))
        % The value is either a single string, or a single numerical value.
        % Do not add {} around it, unless it's one of the exceptional
        % fields that requires that.
        
        field_as_string = string(field_value);
        
        if does_field_need_braces(field_name)
            field_as_string = sprintf('{%s}', field_as_string);
        end
        
        fprintf(fid,'%s = %s\r\n',field_name_for_writing,field_as_string);
    elseif iscell(field_value)
        % The value is a cell array (probably of strings). Add {}.
        
        % We're gonna need it to be one-dimensional
        field_value = field_value(:);
        
        if isempty(field_value)
            % 0 items to write. The standard method will error out in this
            % case, but it's easy to express the result we need, thus:
            fprintf(fid,'%s = {}\r\n',field_name_for_writing);
            
        else
            % We need enough characters to write each string, plus a comma and
            % space after each one, except the last, which we'll just put a }
            % after
            total_characters = sum(strlength(field_value)) + 2*length(field_value) - 1;
            average_characters_per_item = total_characters / length(field_value);
            characters_per_line_max = 80;

            % Ballpark the number of items we can put on each line
            items_per_line = floor(characters_per_line_max / average_characters_per_item);

            value_as_string = '';
            for line_start_index = 1:items_per_line:length(field_value)
                new_line = sprintf('%s, ',field_value{line_start_index:min(line_start_index + items_per_line - 1, length(field_value))});
                value_as_string = sprintf('%s%s\r\n',value_as_string,new_line);
            end

            % replace the comma & space at the end of the string with a } (they're right before the carriage return & newline)
            value_as_string(end-3:end) = ''; % remove last comma & space 
            value_as_string = [value_as_string, '}'];

            if items_per_line >= length(field_value)
                % Everything fit on one line
                % (the ending } is already in value_as_string)
                fprintf(fid,'%s = {%s\r\n',field_name_for_writing,value_as_string);
            else
                % Multiple lines were needed
                fprintf(fid,'%s = {\r\n',field_name_for_writing);
                fprintf(fid,'%s\r\n', value_as_string);
            end
            
        end
        
    else
        fclose(fid); % Cleanin' up
        error('The field %s has a value of unknown type.', field_name);
    end
end 

% Finish and close file
fclose(fid);

end

function fieldname_out = unclean_field_name(fieldname_in)
% Transform a Matlab structure field name into the form it must have to be
% a field in an ENVI header. This undoes
% read_envi_header::clean_field_name, by changing underscores back to
% spaces
fieldname_out = regexprep(fieldname_in,'_',' ');
end

function bool = does_field_need_braces(fieldname)
% Returns true if the field is one of a small number of single-string
% fields that needs to be wrapped in curly braces in the header.

needs_brace_fields = {'description', 'map_info', 'projection_info', 'coordinate_system_string'};
bool = any(strcmp(fieldname, needs_brace_fields));

end

function output_list = permute_strings_array(input_list, desired_order)
% Shuffles around the array of strings _input_list_ so that any of its
% elements that are in _desired_order_ are in the same order as in that
% list and are at the front (with any other elements being at the end in
% some arbitrary order).

output_list = {}; % cell array

% First pass: transfer over any strings that are in desired_order
for i=1:length(desired_order)
    desired_order_element = desired_order{i};
    count_in_input_list = sum(strcmp(desired_order_element, input_list));
    for j=1:count_in_input_list
        output_list{end+1} = desired_order_element;
    end
end

% Second pass: transfer over all strings that aren't in desired_order
for i=1:length(input_list)
    input_list_element = input_list{i};
    is_in_desired_order = any(strcmp(input_list_element, desired_order));
    if ~is_in_desired_order
        output_list{end+1} = input_list_element;
    end
end

% Now output_list should contain all of input_list's elements

end
