function recon = reconstruct_by_cwt(new_dataset, filled_final, wl)
% RECONSTRUCT_BY_CWT_GENERAL Remove etalon noise via CWT masking and ICWT reconstruction
%   Inputs:
%       new_dataset - M x N matrix of input Raman spectra
%       filled_final - binary mask of etalon band (F x N)
%       wl - wavelength axis (1 x N)
%   Output:
%       recon - reconstructed signals (M x N)

    fs = 1 / (wl(3) - wl(2));
    recon = zeros(size(new_dataset));

    for i = 1:size(new_dataset, 1)
        wt = cwt(new_dataset(i,:), 'amor', fs);
        wt_clean = wt .* (1 - filled_final);
        recon(i,:) = smoothdata(icwt(wt_clean, 'SignalMean', mean(new_dataset(i,:))), 'movmedian', 9);
    end
end

