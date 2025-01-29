function [out_matrix, header_data] = read_envi_data( data_filename )
%READ_ENVI_DATA Reads a data file using an accompanying ENVI header file
%   Reads a data file of type 'ENVI Standard','ENVI Spectral Library', or
%   'TIFF' using information found in its ENVI header file. Returns
%   _out_matrix_, the matrix of the file's data, with the dimensions being
%   ordered: lines by samples by bands.
%   Also returns the _header_data_, a structure with the values of all the
%   fields found in the header.

% Get the header data. Header file name is this one's name w/ extension
% replaced
header_filename = get_header_filename(data_filename);
header_data = read_envi_header(header_filename);

data_dims = [header_data.lines,header_data.samples,header_data.bands];

if any(strcmp(header_data.file_type,{'ENVI','ENVI Standard','ENVI Spectral Library'}))
    
    % Test to make sure that the file really contains the amount of data
    % its header claims
    %validate_datafile_size(data_filename, data_dims, header_data.dataType);

    out_matrix = multibandread(data_filename,...
        data_dims,...
        header_data.data_type_matlab,...
        header_data.header_offset,...
        header_data.interleave,...
        header_data.byte_order_matlab);
    
elseif strcmp(header_data.file_type,'TIFF')
    
    out_matrix = imread(data_filename);
    
    % Cast to double, to match the behavior of the ENVI standard part (only
    % ever outputting doubles, regardless of what was in the file)
    out_matrix = double(out_matrix);
    
else
    error('Data file has file_type ''%s'', which this reader doesn''t know how to handle', header_data.file_type);
end

end

function header_filename = get_header_filename(data_filename)
    [pathstr,mainname,ext] = fileparts(data_filename);
    
    if isempty(pathstr)
        pathstr = '.';
    end
    
    hypothetical_header_filename = sprintf('%s.hdr',mainname);
    
    % our guess at the filename might have incorrect case (and this matters
    % in Linux). See if any of the data file's siblings are of the correct
    % form.
    
    directory_contents = getContentsList(pathstr);
    
    files_match_vect = strcmpi(hypothetical_header_filename,directory_contents); % case-insensitive
    
    % error checking
    if sum(files_match_vect)>1
        error('There are multiple candidates for the header file %s in path %s',hypothetical_header_filename,pathstr);
    elseif sum(files_match_vect)==0
        error('There are no candidates for the header file %s in path %s',hypothetical_header_filename,pathstr);
    end
    
    header_filename = fullfile(pathstr,directory_contents{files_match_vect});
    
end

function filenames = getContentsList(pathstr)
    fileStruct = dir(pathstr);
    filenames = cell(length(fileStruct),1);
    for i=1:length(fileStruct)
        filenames{i} = fileStruct(i).name;
    end
end