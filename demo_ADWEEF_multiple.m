% demo_ADWEEF.m - Example script for running ADWEEF algorithm
% This script demonstrates the full denoising pipeline using CWT-based reconstruction
% for dual-wavelength Raman spectroscopy data.

clc; clear; close all;

%% Add function path
addpath('functions');

%% Load wavelength reference files
load('data/wl_671.mat');  % wl_671
load('data/wl_785.mat');  % wl_785
load('data/wl.mat');      % wl

%% Define input data path
load('data/example_multiple_671.mat')
load('data/example_multiple_785.mat')

%% Define parameters
num_components = 5;   % Number of ICA components
loop_num = 5;         % Number of ICA-CWT loops
cutoff_freq = 0.8;    % Cutoff frequency for final Fourier filter

%% Run main pipeline
[recon_671, recon_785] = main_ADWEEF_multiple(new_dataset671, new_dataset785, wl_671, wl_785, wl, num_components, loop_num, cutoff_freq);

%% Plot results for selected depth range
depth_start = 6;  % You may change this depending on your data
disp('Plotting raw spectra');
plot_reconstructed_spectra(new_dataset671, new_dataset785, wl_671, wl_785, depth_start);
disp('Plotting reconstructed spectra');
plot_reconstructed_spectra(recon_671, recon_785, wl_671, wl_785, depth_start);
