function etlon_recon = reconstruct_etalon(wt, filled_final, dataset)
% RECONSTRUCT_ETALON  Reconstruct etalon pattern from masked wavelet coefficients
%
% Inputs:
%   wt           - wavelet transform of signal (F x N)
%   f            - frequency axis corresponding to wavelet transform
%   filled_final - binary mask of etalon frequency band (F x N)
%   dataset      - original signal (M x N), each row is one spectrum
%   loop_times   - current iteration index (for reproducibility if needed)
%
% Output:
%   etlon_recon  - reconstructed etalon signal (1 x N)

    etlon_recon = zeros(1, size(wt, 2));

    etlon_wt = wt .* filled_final;
    etlon_recon = icwt(etlon_wt, 'SignalMean', mean(dataset));
    etlon_recon = etlon_recon - etlon_recon(1);
end
