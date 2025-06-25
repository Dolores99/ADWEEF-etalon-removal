function recon = reconstruct_by_subtraction(new_dataset, etalon_recon, wl)
% RECONSTRUCT_BY_SUBTRACTION Remove etalon noise via direct subtraction/addition
%   Inputs:
%       new_dataset   - M x N matrix of input Raman spectra
%       etalon_recon  - 1 x N vector of reconstructed etalon noise
%       wl            - wavelength axis (1 x N)
%   Output:
%       recon - reconstructed signals (M x N)

    recon = zeros(size(new_dataset));
    fs = 1 / (wl(3) - wl(2));
    [~, para] = max(wl .* (wl <= 870));

    for i = 1:size(new_dataset, 1)
        max_val = max(new_dataset(i, :));
        min_val = min(new_dataset(i, :));
        norm_signal = (new_dataset(i, :) - min_val) / (max_val - min_val);

        if std(norm_signal - etalon_recon) < std(norm_signal + etalon_recon)
            cleaned = norm_signal - etalon_recon;
        else
            cleaned = norm_signal + etalon_recon;
        end

        cleaned(para:end) = smoothdata(cleaned(para:end), 'movmedian', 9);
        recon(i,:) = cleaned * (max_val - min_val) + min_val;
    end
end
