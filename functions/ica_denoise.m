function [unmixed, weights, noise_indices] = ica_denoise(dataset, num_components)
%ICA_DENOISE Perform ICA decomposition and detect noise components
%
%   [unmixed, weights, noise_indices] = ica_denoise(dataset, num_components)
%
%   Input:
%       dataset         - (M x N) Raman data, each row is one spectrum
%       num_components  - number of ICA components
%
%   Output:
%       unmixed         - (N x num_components) ICA components (normalized)
%       weights         - (M x num_components) ICA transform weights
%       noise_indices   - indices of components considered as noise

    % Transpose to (N x M) before applying RICA
    dataset_T = dataset';

    % Perform ICA via RICA
    Mdl = rica(dataset_T, num_components, 'IterationLimit', 1e4);
    weights = Mdl.TransformWeights;
    unmixed = transform(Mdl, dataset_T);

    % Detect noise components (weights are all-positive or all-negative)
    noise_indices = [];
    for col = 1:size(unmixed, 2)
        col_weights = weights(:, col);
        if all(col_weights >= 0) || all(col_weights <= 0)
            noise_indices(end+1) = col;
        end
    end

    % Normalize each component to [0, 1]
    for col = 1:size(unmixed, 2)
        col_data = unmixed(:, col);
        unmixed(:, col) = (col_data - min(col_data)) / (max(col_data) - min(col_data));
    end
end
