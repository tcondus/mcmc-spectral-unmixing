function [data_filename, header_filename] = write_envi_data(data,out_filename_requested,interleave,varargin)
%WRITE_ENVI_DATA writes file with data and a corresponding ENVI header
% Takes a matrix of data, a filename & interleave for the data, and maybe
% additional arguments to be passed to save_envi_header (giving information
% like description text or wavelengths to be saved in the header). See that
% function's notes for a listing of understood options.
% RETURNS: the name it in fact used for the datafile (which might be
% different from the argument)

% Possible parameters in varargin:
% - 'inherit_from': a structure of ENVI header fields, as produced by
% read_envi_header. The values in this structure, if present, will be taken
% as the defaults. They can be overridden by any of the parameters listed
% below. This argument can also be used to add arbitrary fields to the
% header
% - 'desc': text to use, if not there then generic msg
% - 'type': the type to translate into ENVI format. If absent, than use
% the type of data param
% - 'wavelengths': a vector to use. If absent, omit section altogether.
% - 'default bands': an RGB vector to use (length 3). If absent, omit.
% - 'mapProjectionInfo': an object describing map projection
% information relavent to the saved data. If absent, omit.
% - 'badbandslist': a vector to put as the bad bands list (bbl). If
% absent, omit section.
% - 'band_names': a list of names for the bands given as a cell array.
% Maps into the header record type labeled "band names".

[pathname, basename, ext] = fileparts(out_filename_requested);
% If there's no extension, append a generic one based on interleave
% mode. Otherwise, leave the filename alone
if strcmp(ext, '')
    data_filename = strcat(out_filename_requested, '.', interleave);
    header_filename = strcat(out_filename_requested, '.hdr');
else
    data_filename = out_filename_requested;
    header_filename = fullfile(pathname,strcat(basename, '.hdr'));
end

args = parseParameters(varargin);

header_data = [];

if isfield(args, 'inherit_from')
    % There's a struct we should use the properties from (with other
    % command line arguments overriding its values)
    header_data = args.inherit_from;
end

% Use the command line arguments given to this function to override
% values in the inherited header info
if isfield(args, 'desc')
    header_data.description = args.desc;
end
if isfield(args, 'type')
    header_data.data_type_matlab = args.type;
    % TODO: this one's a bit irregular! We shouldn't pretend that the
    % data_type thing's gonna be there, or be right, after this
    % assignment for sure.
end
if isfield(args, 'wavelengths')
    header_data.wavelength = args.wavelengths;
    if ~isfield(header_data, 'wavelength_units')
        % If they included wavelengths but didn't set units, we're going to
        % assume it's nanometers (a previous version of the header-writing
        % code assumed this)
        header_data.wavelength_units = 'Nanometers';
    end
end
if isfield(args, 'default_bands')
    switch length(args.default_bands)
        case 3
            % Good default bands
            header_data.default_bands = args.default_bands;
        case 0
            % There's nothing there. Quietly don't bother with it.
            % Use value in the inherit_from stuff, if there is one.
        otherwise
            error('Vector of default bands has length %i (should be 3, with entries being RGB respectively)', length(args.default_bands));
    end
end
if isfield(args, 'mapProjectionInfo')
    if ~isempty(args.mapProjectionInfo)
        header_data.map_info = args.mapProjectionInfo.mapInfo;
        header_data.projection_info = args.mapProjectionInfo.projInfo;
        header_data.coordinate_system_string = args.mapProjectionInfo.coSys;
    end
end
if isfield(args, 'badbandslist')
    header_data.bbl = args.badbandslist;
end
if isfield(args, 'band_names')
    header_data.band_names = args.band_names;
end

% Also, we don't care a wit about the interleaving in the inherit_from
% stuff. (Backwards compatibility is a helluva drug. I'd never build
% it like this if I were starting now.)
header_data.interleave = interleave;

% Also, any sizes given in the inherit_from data are unimportant. We
% will just use the actual data's size
[lines,samples,bands] = size(data);
header_data.samples = samples;
header_data.lines = lines;
header_data.bands = bands;

% Also, by fiat we're forcing header offset and byte order to a
% standard orientation.
header_data.header_offset = 0; % bytes
header_data.byte_order_matlab = 'ieee-le';
% TODO: this won't set the numeric byte order! That's gotta be done
% somewhere before we save!

% Assuming that this is ENVI standard (for now)
% TODO
header_data.file_type = 'ENVI Standard';

% write data
if ~(isfield(header_data, 'data_type_matlab') || isfield(header_data, 'data_type'))
    % no requested precision: use the storage class of the data in
    % Matlab
    header_data.data_type_matlab = class(data);
end

if ~any(strcmp(header_data.interleave,{'bsq','bil','bip'}))
    error('error in write_envi_data: invalid interleaving method ''%s''',header_data.interleave);
end

% Actually save the data
if isempty(data)
    % No data: simply create empty file and emit warning
    f = fopen(data_filename, 'w');
    fclose(f);
    fprintf('WARNING: writing ENVI data file with size 0 at %s\n', data_filename);
else
    multibandwrite(data, data_filename, interleave, ...
        'precision', header_data.data_type_matlab, ...
        'offset', header_data.header_offset, ...
        'machfmt', header_data.byte_order_matlab);
end

% write header
write_envi_header(header_filename, header_data);

end
