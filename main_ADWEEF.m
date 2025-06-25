function [recon_671, recon_785] = main_ADWEEF(csv_path, wl_671, wl_785, num_components, loop_num, cutoff_freq)
% MAIN_ADWEEF_RECONSTRUCTION Complete automated pipeline for Etalon removal
%
% Inputs:
%   csv_path       - full path to .csv Raman data file
%   wl_671         - 1 x N wavelength axis for 671 nm
%   wl_785         - 1 x N wavelength axis for 785 nm
%   num_components - number of ICA components (e.g. 5)
%   loop_num       - number of ICA-CWT iterations (e.g. 5)
%   cutoff_freq    - frequency cutoff for final Fourier filter (e.g. 0.8)
%
% Outputs:
%   recon_671      - final denoised signals for 671 nm (M x N)
%   recon_785      - final denoised signals for 785 nm (M x N)

    %% Step 1: Load Data
    [wl, dataset671, dataset785] = load_raman_data(csv_path);
    new_dataset671 = dataset671;
    new_dataset785 = dataset785;

    %% Iterative ICA-CWT Reconstruction
    for loop_times = 1:loop_num
        % ICA Denoising
        [unmixed_671, weights_671, noise_idx_671] = ica_denoise(new_dataset671, num_components);
        [unmixed_785, weights_785, noise_idx_785] = ica_denoise(new_dataset785, num_components);

        % CWT and Mask Generation
        [wt_671, f_671] = cwt(unmixed_671(:, noise_idx_671(1)), 'amor', 1 / (wl(3)-wl(2)));
        [wt_785, f_785] = cwt(unmixed_785(:, noise_idx_785(1)), 'amor', 1 / (wl(3)-wl(2)));
        filled_final = build_mask_from_cross_cwt(wt_671, wt_785, f_671, wl);

        % Reconstruct signals by masking wavelet coefficients and inverse CWT
        recon_671 = reconstruct_by_cwt(new_dataset671, filled_final, wl);
        recon_785 = reconstruct_by_cwt(new_dataset785, filled_final, wl);

        new_dataset671 = recon_671;
        new_dataset785 = recon_785;
    end

    %% Final low-pass Fourier filtering
    recon_671 = fourier_filter(recon_671, wl, cutoff_freq);
    recon_785 = fourier_filter(recon_785, wl, cutoff_freq);
end
