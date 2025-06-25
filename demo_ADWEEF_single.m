% DEMO_ADWEEF_SINGLE  Run single-depth ADWEEF subtraction pipeline

% Add path to functions if needed
% addpath('functions');

% Load wavelength vectors
load('data/wl_671.mat');
load('data/wl_785.mat');
load('data/wl.mat');

% CSV path
csv_path = 'data/example_single.csv';

% Loop count
loop_num = 5;

% Run single-depth Etalon removal
[recon_671, recon_785] = main_ADWEEF_single(csv_path, wl_671, wl_785, wl, loop_num);

% Plot results
figure;
subplot(2,1,1);
plot(wl_671, recon_671, 'r', 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Intensity');
title('Reconstructed Signal - 671 nm');

subplot(2,1,2);
plot(wl_785, recon_785, 'b', 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Intensity');
title('Reconstructed Signal - 785 nm');
