%% run_adaptive_doppler_estimation.m
% -------------------------------------------------------------------------
% References:
% [1] S.M.M. Tabatabaei Majd, L. Eslami, B. Mohammadzadeh Asl,
%     "Adaptive Doppler blood flow estimation in ultrasound with enhanced spectral resolution and contrast using limited observation windows,"
%     Ultrasonics, vol. 154, 2025, 107678.
%     https://doi.org/10.1016/j.ultras.2025.107678
%
% [2] S.M.M. Tabatabaei Majd, B. Mohammadzadeh Asl,
%     "Adaptive Spectral Doppler Estimation Based on the Modified Amplitude Spectrum Capon,"
%     IEEE Trans. UFFC, vol. 68, no. 5, pp. 1664–1672, May 2021.
%     https://doi.org/10.1109/TUFFC.2020.3044774
% =========================================================================
% Adaptive Doppler Spectral Estimation on Femoral and Laminar Flow Data
%
% Implements Welch, Capon, Projection-Capon, MASC, and HQASC estimators.
% Authors:
% Seyed Mohammad Mahdi Tabatabaei Majd, Leila Eslami
% =========================================================================

clear; clc; close all;

%% --- Step 1: Load Doppler IQ Data ---
[filename, pathname] = uigetfile('*.mat', 'Select Doppler dataset');
if isequal(filename, 0)
    error('No file selected. Exiting...');
end
load(fullfile(pathname, filename));  % expects 'data' variable
fprintf('Loaded data from %s\n', filename);

%% --- Step 2: Define Physical and Acquisition Parameters ---
f0 = 5e6;
fs = 20e6;
fprf = 15e3;
angle_deg = 60;
c = 1540;
dB_range = 40;

%% --- Step 3: Analysis Setup ---
OWs = [2 128];
Averaging_Kernels = [4 256];
M = 500;
dN = 1;
phi = f0 / fs;
w_t = 2 * pi * phi;
[depths, frames] = size(data);
allSpectra = cell(1, length(OWs));

%% --- Step 4: Select Folder for Output ---
savePath = uigetdir(pwd, 'Select folder to save output results');
if savePath == 0
    error('No folder selected. Exiting...');
end
processed_name = erase(filename, '.mat');

%% --- Step 5: Main Spectral Estimation Loop ---
for idx = 1:length(OWs)
    N = OWs(idx);
    L = Averaging_Kernels(idx);
    nSpec = floor((frames - L) / dN) + 1;

    if N == 128
        methods = {'Welch-Rect', 'Welch-Hamming'};
        Spectra = zeros(M, nSpec, 2);
    else
        methods = {'Welch', 'Capon', 'Pr.Capon', 'MASC', 'HQASC'};
        Spectra = zeros(M, nSpec, 5);
    end

    for i = 0:nSpec - 1
        fprintf('OW = %d | Spectrum %d of %d\n', N, i+1, nSpec);
        segment = data(:, i*dN + 1 : i*dN + L);

        if N == 128
            [P_rect, P_ham] = REF_Spectral_Ham(segment, N, M);
            Spectra(:, i+1, 1) = P_rect;
            Spectra(:, i+1, 2) = P_ham;
        else
            [P1, P2, P3, P4, P5] = adaptive_spectral_estimators(segment, N, M, w_t);
            Spectra(:, i+1, :) = cat(3, P1, P2, P3, P4, P5);
        end
    end

    allSpectra{1, idx} = Spectra;
end

%% --- Step 6: Save Spectra Output ---
outFile = fullfile(savePath, ['allSpectra_' processed_name '.mat']);
save(outFile, 'allSpectra');
fprintf('Saved spectra to %s\n', outFile);

%% --- Step 7: Plot Spectrograms ---
doPlot = true;
if doPlot
    flipData = false;
    setup = struct('f0', f0, 'fs', fs, 'fprf', fprf, 'c', c, 'theta', angle_deg);

    for t = 1:length(OWs)
        N = OWs(t);
        Spec = allSpectra{1, t};
        T = size(Spec, 2);
        timeAxis = (0:T - 1) * dN / fprf;
        freqAxis = fprf / M * (-(M/2):(M/2 - 1));
        velAxis = c / 2 * freqAxis / f0 / cosd(angle_deg);

        for k = 1:size(Spec, 3)
            titleStr = sprintf('Method: %s | OW: %d', methods{k}, N);
            extra.figTitle = titleStr;
            plotSpectrogram(Spec(:, :, k), timeAxis, velAxis, dB_range, flipData, extra);
        end
    end
end
