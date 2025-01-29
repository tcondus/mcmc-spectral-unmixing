function [ headerData ] = read_envi_header( filename )
%READ_ENVI_HEADER Extracts items from an ENVI header
%   Extracts all information from an ENVI header. If the function succeeds,
%   and the data file is a standard flat array, this is guaranteed to
%   include 'samples', 'lines', 'bands', 'interleave', 'header_offset',
%   'byte_order', 'data_type', 'file_type'. If the data file is a TIFF, the
%   output is guaranteed to have all those fields except 'interleave',
%   'header_offset', and 'byte_order'.
%
%   All field names in the
%   results here are the same as in the header itself, with spaces replaced
%   by underscores.
%
%   Note that we add some derived fields that aren't present in the file:
%    - data_type_matlab: a string based on data_type that is ready to be
%    given to multibandread
%    - byte_order_matlab: a string based on byte_order that is ready to be
%    given to multibandread (except for TIFFs)

fid = fopen(filename);

if fid==-1
    error('The file %s does not exist!\n',filename);
end

allLines = {};

% read lines until file is consumed
errnum = 0;
while errnum == 0
    nextLine = fgetl(fid);
    
    allLines{end+1} = remove_comments(nextLine);
    
    [~,errnum] = ferror(fid); % -4 is EOF reached
end
fclose(fid);

num_lines = length(allLines);

headerData = [];
waiting_for_group_end = false;
current_group_name = '';
for i=1:num_lines
    thisLine = allLines(i);
    % There are 3 things that can be in a line of an ENVI header, once
    % we've removed the comments:
    %
    % 1) The start of a Key-value pair (which might have a simple value or
    % a group in curly braces. If it's a curly-brace group, it may or may
    % not end on this line.)
    % 2) Values in the middle/end of a group set off by curly braces
    % 3) Something semantically void, like a line of nothing (or perhaps
    % spaces), or the characters 'ENVI' found in the first line
    %
    % Anything else is an error. (Although, at this time, we don't bother
    % protecting
    
    if waiting_for_group_end
        % We are in the midst of a group. We will probably be adding stuff
        % to the current group.
        
        % Add the text (except the line break) into the field's text in our
        % structure
        
        thisLineText = thisLine{1};
        
        fieldsofar = headerData.(current_group_name);
        headerData.(current_group_name) = strcat(fieldsofar, thisLineText);
        
        if contains(thisLineText, '}')
            waiting_for_group_end = false;
            
            % Get rid of the trailing }, which there surely is
            data_str = headerData.(current_group_name);
            data_str = strip_trailing_curlybrace(data_str);
            
            % Break the data down into an array of appropriate type
            headerData.(current_group_name) = group_to_array(data_str, current_group_name);
        end
        
    else
        % We are not in a group. Look for a simple value, or the start of a
        % group.
        if contains(thisLine, '=')
            if contains(thisLine, '{')
                % A group start, and perhaps end
                captured = regexp(thisLine, '(.*?)\W*=\W*\{(.*)', 'tokens');
                captured = captured{1}{1}; % Shed some outer layers of the regexp result structure
                fieldname = captured{1};
                fieldfirstvalues = captured{2}; % Might have a } at the end
                
                current_group_name = clean_field_name(fieldname); % Now suitable for saving in Matlab
                
                if contains(thisLine, '}')
                    % The group ends on this line, too.
                    
                    % Get rid of the trailing }
                    fieldfirstvalues = strip_trailing_curlybrace(fieldfirstvalues);
                    
                    % Break the data down into an array of appropriate type
                    headerData.(current_group_name) = group_to_array(fieldfirstvalues, current_group_name);
                else
                    % There will be more lines to the group
                    
                    % Start an entry for this field in our structure
                    headerData.(current_group_name) = fieldfirstvalues;
                    
                    waiting_for_group_end = true;
                end
                
            else
                % A simple value
                captured = regexp(thisLine, '(.*?)\s*=\s*(.*)', 'tokens');
                captured = captured{1}{1}; % Shed some outer layers of the regexp result structure
                fieldname = captured{1};
                fieldvalue_string = strip(captured{2}); % Forget leading & trailing spaces of the value
                
                current_group_name = clean_field_name(fieldname); % Now suitable for saving in Matlab
                
                fieldvalue_number = str2double(fieldvalue_string);
                
                if isnan(fieldvalue_number)
                    % Not a number. We will save the string, then
                    headerData.(current_group_name) = fieldvalue_string;
                else
                    % Since it is a number, save it as a number
                    headerData.(current_group_name) = fieldvalue_number;
                end
                
            end
        else
            % There's not a key-value pair here. We'll just assume it's not
            % important, then.
            
            % (Intentionally doing nothing)
        end
    end
    
end

if ~all(isfield(headerData, {'samples', 'lines', 'bands', 'data_type', 'file_type'}))
    error('Header lacks at least one of these essential fields: samples, lines, bands, data type, file type.');
end

% get name of data type from the numeric code given in the header
switch headerData.data_type
    case 1
        headerData.data_type_matlab = 'uint8';
    case 2
        headerData.data_type_matlab = 'int16';
    case 3
        headerData.data_type_matlab = 'int32';
    case 4
        headerData.data_type_matlab = 'float'; % not 'single', because these types are meant for multibandread
    case 5
        headerData.data_type_matlab = 'double';
    case 12
        headerData.data_type_matlab = 'uint16';
    case 13
        headerData.data_type_matlab = 'uint32';
    case 14
        headerData.data_type_matlab = 'int64';
    case 15
        headerData.data_type_matlab = 'uint64';
    otherwise
        error('Header has unknown data type %i. You might consider modifying this function to support this type.', headerData.data_type);
end

if any(strcmp(headerData.file_type,{'ENVI','ENVI Standard','ENVI Spectral Library'}))
    % This is one of the standard ENVI flat files.
    
    % Check to make sure that the essential elements are present in the header
    if ~all(isfield(headerData, {'interleave', 'header_offset', 'byte_order'}))
        error('Header for a standard ENVI data file lacks at least one of these essential fields: interleave, header offset, byte order.');
    end

    switch headerData.byte_order
        case 0
            headerData.byte_order_matlab = 'ieee-le';
        case 1
            headerData.byte_order_matlab = 'ieee-be';
        otherwise
            error('Header has unknown byte order %i. Only 0 (little-endian) and 1 (big-endian) are valid.', headerData.byte_order);
    end

elseif strcmp(headerData.file_type,'TIFF')
    % TIFFs give us slightly different data needs:
    % - we don't need interleave, header offset, or byte order
    
else
    error('Header is for file type ''%s'', which this reader doesn''t know how to handle', headerData.file_type);
end

end

function cleaned_line = remove_comments(text_line)
% Removes ENVI comments, which are lines in which the first non-whitespace
% character is a semicolon. (Semicolons elsewhere don't cause the text
% after them to be a comment.)
% Requires that there be no trailing newline in the input _text_line_.
cleaned_line = regexprep(text_line,'^\s*;.*$','');
end

function fieldname_out = clean_field_name(fieldname_in)
% Transform a header field name so that it can be a field in a Matlab
% structure (i.e., no spaces)
fieldname_out = regexprep(fieldname_in,' ','_');
end

function out_str = strip_trailing_curlybrace(in_str)
% Given a string that has some values at the end of a group, and then a
% closing curly brace, returns the string up to (but not including!)
% that closing brace.

% The .*? is lazy, so that we only get things up to the first } we
% found. The case where there could be multiple will requires some more
% thought.
captured = regexp(in_str,'(.*?)\}','tokens');
out_str = captured{1}{1};
end

function out_array = group_to_array(group_string, field_name)
% Given a string that is the comma-separated contents of a group in an
% ENVI label (without the enclosing curly braces), produces an
% appropriate array. If all the values in the group are numeric, that
% array is an array of those numbers as doubles. If they're not, it's a
% cell array of strings.

% NOTE: if the field_name matches one of a select number of known fields,
% it will not be split! These are fields whose values are special strings
% that are not actually meant to be taken as a group.

nosplit_fields = {'description', 'map_info', 'projection_info', 'coordinate_system_string'};

if any(strcmp(field_name, nosplit_fields))
    % No splitting? Just strip, then.
    out_array = strip(group_string);
    return;
end

% Split along commas
contents_strings_raw = split(group_string, ',');

num_values = length(contents_strings_raw);
contents_strings_stripped = cell(1, num_values);
contents_numeric = zeros(1, num_values);
are_all_numeric = true;

for i= 1:num_values
    % Strip each value, and try to convert them into numbers
    
    stripped = strip(contents_strings_raw{i});
    contents_strings_stripped{i} = stripped;
    
    contents_numeric(i) = str2double(stripped);
    
    if isnan(contents_numeric(i))
        are_all_numeric = false;
    end
end

% Are all of the values numbers? If so, we will return that.
if are_all_numeric
    out_array = contents_numeric;
else
    out_array = contents_strings_stripped;
end
end
