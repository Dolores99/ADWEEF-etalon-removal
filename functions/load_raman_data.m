function [wl, dataset671, dataset785] = load_raman_data(filepath)
% LOAD_RAMAN_DATA Load dual-wavelength Raman data from CSV
%
%   [wl, dataset671, dataset785] = load_raman_data(filepath)
% 
%   Input:
%       filepath - full path to the CSV file
%
%   Output:
%       wl         - wavelength axis (1 x N)
%       dataset671 - Raman signals under 671 nm (M x N)
%       dataset785 - Raman signals under 785 nm (M x N)

    data = readmatrix(filepath);
    wl = data(1, :);                    % First row: wavelength axis
    dataset671 = data(4:2:end, :);      % Even-indexed rows (2, 4, ...) as 671 nm data
    dataset785 = data(5:2:end, :);      % Odd-indexed rows (3, 5, ...) as 785 nm data
end
