function [recon_671, recon_785] = main_ADWEEF_single(csv_path, wl_671, wl_785, wl_all, loop_num, cutoff_freq)
% MAIN_SUBTRACTION_SINGLE  Etalon removal by direct subtraction (no ICA), single-depth spectra
%
% Inputs:
%   csv_path  - path to 5-line CSV (wl, wl, wl, 671 signal, 785 signal)
%   wl_671    - wavelength axis for 671 nm
%   wl_785    - wavelength axis for 785 nm
%   wl_all    - global wavelength axis
%   loop_num  - number of refinement loops (e.g., 5)
%
% Output:
%   recon_671 - reconstructed signal (1 x N)
%   recon_785 - reconstructed signal (1 x N)

    %% Load single-depth data
    data = readmatrix(csv_path);
    wl = data(1, :);
    signal_671 = data(4, :);
    signal_785 = data(5, :);

    % For consistency with multi-depth code, wrap into (1 x N)
    dataset671 = signal_671;
    dataset785 = signal_785;

    %% Iterative subtraction
    for loop_times = 1:loop_num
        % Continuous Wavelet Transform
        fs = 1 / (wl(3) - wl(2));
        [wt_671, f_671] = cwt(dataset671, 'amor', fs);
        [wt_785, f_785] = cwt(dataset785, 'amor', fs);

        % Build mask via CWT cross-spectrum
        filled_final = build_mask_from_cross_cwt(wt_671, wt_785, f_671, wl);

        % Reconstruct via subtraction (add or subtract depending on std)
        recon_671 = reconstruct_by_subtraction(dataset671, filled_final, wl);
        recon_785 = reconstruct_by_subtraction(dataset785, filled_final, wl);

        dataset671 = recon_671;
        dataset785 = recon_785;
    end

    %% Final low-pass Fourier filtering
    recon_671 = fourier_filter(recon_671, wl, cutoff_freq);
    recon_785 = fourier_filter(recon_785, wl, cutoff_freq);
end
