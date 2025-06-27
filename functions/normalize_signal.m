function norm_data = normalize_signal(data)
% NORMALIZE_SIGNAL Normalize each signal to [0, 1]
%
%   norm_data = normalize_signal(data)
%
%   Input:
%       data      - (M x N) matrix (each row a signal), or 1D row/column vector
%
%   Output:
%       norm_data - same size as input, values normalized to [0, 1]

    sz = size(data);
    if isvector(data)
        % If 1D vector (row or column), normalize it directly
        min_val = min(data(:));
        max_val = max(data(:));
        if max_val ~= min_val
            norm_data = (data - min_val) / (max_val - min_val);
        else
            norm_data = zeros(size(data));
        end
    else
        % If 2D matrix, normalize each row to [0, 1]
        norm_data = zeros(sz);
        for i = 1:sz(1)
            row = data(i, :);
            min_val = min(row);
            max_val = max(row);
            if max_val ~= min_val
                norm_data(i, :) = (row - min_val) / (max_val - min_val);
            else
                norm_data(i, :) = zeros(1, sz(2));
            end
        end
    end
end
