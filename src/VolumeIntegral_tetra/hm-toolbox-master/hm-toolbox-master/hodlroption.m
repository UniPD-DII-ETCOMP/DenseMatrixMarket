function opt = hodlroption(key, value)
%HODLROPTION Set or get an option for the hodlr toolbox.
%
% Valid options are:
%   'block-size': Integer representing the minimum block size.
%   'threshold': Value used for off-diagonal truncation.
%   'compression': Can be either 'qr' or 'svd', and selects the
%       method used to compress unstructured dense matrices when calling
%       hodlrA = hodlr(A);
%   'norm': The norm used for truncation. Can be either 2, or 'fro', and
%       the truncation is performed relative to the norm of the entire
%       matrix. 
%
% The special option HODLROPTION('clear') can be used to load all the
% default values, clearing all the previous HODLROPTION commands. 

global hodlr_block_size
global hodlr_threshold
global hodlr_compression
global hodlr_norm

if strcmp(key, 'clear')
    if exist('value', 'var')
        error('Specifying a value is unsupported for the special option "clear"');
    end
    
    clear hodlr_block_size;
    clear hodlr_threshold;
    clear hodlr_compression;
    clear hodlr_norm;
    
    return;
end

if isempty(hodlr_compression)
	hodlr_compression = 'qr';
end

if isempty(hodlr_block_size)
    hodlr_block_size = 256;
end

if isempty(hodlr_threshold)
    hodlr_threshold = 1e-12;
end

if isempty(hodlr_norm)
    hodlr_norm = 2;
end

if ~exist('key', 'var')
    error('Please specify a key');
end

if ~exist('value', 'var')
    switch key
        case 'block-size'
            opt = hodlr_block_size;
        case 'threshold'
            opt = hodlr_threshold;
		case 'compression'
			opt = hodlr_compression;
        case 'norm'
            opt = hodlr_norm;
        otherwise
            error('Unsupported option specified');
    end
else
    switch key
        case 'block-size'
            if value <= 2
                error('minimum block size must be at least 3');
            else
                hodlr_block_size = value;
            end
        case 'threshold'
            if value < 0
                error('threshold has to be positive');
            else
                hodlr_threshold = max(eps, value);
			end
		case 'compression'
			if ~strcmp(value, 'qr') && ~strcmp(value, 'svd') && ...
					~strcmp(value, 'lanczos')
				error('Invalid value for dense-compression');
			else
				hodlr_compression = value;
            end
        case 'norm'
            if ( ischar(value) && ~strcmp(value, 'fro') ) && value ~= 2
                error('Invalid valud for norm');
            else
                hodlr_norm = value;
            end
        otherwise
            error('Unsupported option specified');
    end
    
end

