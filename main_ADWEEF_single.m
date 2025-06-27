function [recon_671, recon_785] = main_ADWEEF_single(signal_671, signal_785, wl, loop_num, cutoff_freq)
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

    
    %% For consistency with multi-depth code, wrap into (1 x N)
    dataset671 = normalize_signal(signal_671');
    dataset785 = normalize_signal(signal_785');

    %% Iterative subtraction
    for loop_times = 1:loop_num
        % Continuous Wavelet Transform
        fs = 1 / (wl(3) - wl(2));
        [wt_671, f_671] = cwt(dataset671, 'amor', fs);
        [wt_785, f_785] = cwt(dataset785, 'amor', fs);

        % Build mask via CWT cross-spectrum
        filled_final = detect_etalon_band(wt_671, f_671, wt_785, f_785, wl, f_671, loop_times);
        if isempty(filled_final)
            break;  % 提前结束循环
        end

        % Reconstruct Etalon
        etlon_recon_671(loop_times,:) = reconstruct_etalon(wt_671, filled_final, dataset671);
        etlon_recon_785(loop_times,:) = reconstruct_etalon(wt_785, filled_final, dataset785);

        % Reconstruct via subtraction (add or subtract depending on std)
        recon_671 = reconstruct_by_subtraction(dataset671', etlon_recon_671(loop_times,:), wl);
        recon_785 = reconstruct_by_subtraction(dataset785', etlon_recon_785(loop_times,:), wl);

        dataset671 = recon_671';
        dataset785 = recon_785';
    end

    %% Final low-pass Fourier filtering
    recon_671 = fourier_filter(recon_671, wl, cutoff_freq);
    recon_785 = fourier_filter(recon_785, wl, cutoff_freq);
end
