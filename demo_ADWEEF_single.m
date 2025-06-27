% DEMO_ADWEEF_SINGLE  Run single-depth ADWEEF subtraction pipeline
clc;close all; clearvars;
% Add path to functions if needed
addpath('functions');

% Load wavelength vectors
load('data/wl_671.mat');
load('data/wl_785.mat');
load('data/wl.mat');

% Load single-depth data
signal_671 = load('data/example_single_671.mat').dataset_671;
signal_785 = load('data/example_single_785.mat').dataset_785;

% Loop count
loop_num = 5;

%
cutoff_freq = 0.8;

% Run single-depth Etalon removal
[recon_671, recon_785] = main_ADWEEF_single(signal_671, signal_785, wl, loop_num, cutoff_freq);

% Plot results
figure;
subplot(2,1,1);
plot(wl_671, signal_671, 'r', 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Intensity');
xlim([2800 3800]);
title('Raw Signal - 671 nm');

subplot(2,1,2);
plot(wl_785, signal_785, 'b', 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Intensity');
xlim([400 1800]);
title('Raw Signal - 785 nm');

figure;
subplot(2,1,1);
plot(wl_671, recon_671, 'r', 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Intensity');
xlim([2800 3800]);
title('Reconstructed Signal - 671 nm');

subplot(2,1,2);
plot(wl_785, recon_785, 'b', 'LineWidth', 1.5);
xlabel('Wavelength (nm)'); ylabel('Intensity');
xlim([400 1800]);
title('Reconstructed Signal - 785 nm');
