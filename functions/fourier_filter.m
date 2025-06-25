function filtered = fourier_filter(signal, wl, cutoff_max)
% FOURIER_FILTER Apply low-pass Fourier filter to suppress high-frequency noise
%
%   filtered = fourier_filter(signal, wl, cutoff_max)
%
%   Inputs:
%       signal      - M x N matrix of Raman spectra
%       wl          - 1 x N wavelength axis
%       cutoff_max  - maximum frequency to preserve (e.g., 0.8)
%
%   Output:
%       filtered    - filtered Raman spectra (M x N)

    [M, N] = size(signal);
    fs = 1 / (wl(3) - wl(2));
    freq = fs * (0:(N/2)) / N;
    filtered = zeros(size(signal));

    for i = 1:M
        sig_fft = fft(signal(i, :));
        filter_mask = ones(1, N);
        filter_mask(freq > cutoff_max) = 0;
        filter_mask(N:-1:N/2+2) = filter_mask(2:N/2); % Mirror for real signal

        sig_fft_filtered = sig_fft .* filter_mask;
        filtered(i, :) = ifft(sig_fft_filtered, 'symmetric');
    end
end
