% test_velocity_estimation.m
% -------------------------------------------------------------------------
% Monte Carlo simulation for Doppler blood velocity estimation using
% multiple spectral methods under various SNR levels.
%
% Methods evaluated: Welch, Capon, Pr.Capon, MASC, HQASC
% Authors: Seyed Mohammad Mahdi Tabatabaei Majd, Leila Eslami
% -------------------------------------------------------------------------

clc; close all; clear;
disp(['Date and time: ' datestr(now)])

%% --- Simulation Parameters ---
whatTest = 1;
SNR = [-20 -15 -10 -5 0 5];      % SNR values in dB
MC = 200;                       % Monte Carlo iterations
saveData = false;
doPlot = true;

P = 500;                        % Frequency points per spectrum

rng(SNR(1) * MC);               % Set reproducible seed

%% --- Physical Constants ---
f_s = 20e6;                     % Sampling frequency [Hz]
f_c = 5e6;                      % Center frequency [Hz]
w_c = 2*pi*f_c;
phi = w_c / f_s;
c = 1540;                      % Speed of sound [m/s]
f_prf = 15e3;                   % Pulse repetition frequency [Hz]
v_z = 0.5;                      % True axial velocity [m/s]
psi = -2 * w_c / (c * f_prf) * v_z;
velTrue = v_z;

%% --- Simulation Settings ---
K = 33;                         % Number of fast-time samples
Ns = 20;                        % Number of slow-time samples (flow shots)
nbrOfBlocks = 1;
no_filt = Ns / 2;
flowShots = ones(1, Ns);
flowShotsAll = repmat(flowShots, 1, nbrOfBlocks);
NAll = length(flowShotsAll);
NVecAll = 0:NAll-1;
freqVecShifted = ((0:P-1) - P/2) / P;
velVec = 2 * pi * freqVecShifted * (c * f_prf) / (-2 * w_c);

algNames = {'Welch', 'Capon', 'Pr.Capon', 'MASC', 'HQASC'};
nbrOfAlgs = length(algNames);
saveStr = sprintf('MSE_vs_SNR_MC%d_Ns_%d', MC, Ns);
disp(['Save string: ' saveStr])

%% --- Generate Noise-Free Signal ---
YNoNoise = zeros(K, NAll);
z = exp(1i * psi * NVecAll);
for k = 1:K
    YNoNoise(k,:) = exp(1i * phi * (k - 1)) * z;
end
sigPower = mean(abs(YNoNoise(:)).^2);

%% --- Determine Noise Variance ---
snrLen = length(SNR);
sigma2 = sigPower * 10.^(-SNR/10);

%% --- Estimation Loop ---
ampEst = zeros(P, MC, snrLen, nbrOfAlgs);
velEst = zeros(MC, snrLen, nbrOfAlgs);
startTimeTot = cputime;

for k0 = 1:snrLen
    disp(['SNR: ' num2str(SNR(k0)) ' dB'])

    for mc = 1:MC
        if ~mod(mc, 10)
            fprintf('MC iteration %d / %d\n', mc, MC)
        end

        % Generate noise and apply
        noise = sqrt(sigma2(k0)/2) * (randn(K, NAll) + 1i * randn(K, NAll));
        Y = YNoNoise + noise;

        % Spectral estimation
        [ampEst(:,mc,k0,1), ampEst(:,mc,k0,2), ampEst(:,mc,k0,3),
         ampEst(:,mc,k0,4), ampEst(:,mc,k0,5)] = adaptive_spectral_estimators(Y, no_filt, P, phi);

        % Peak-based velocity estimation
        for ii = 1:nbrOfAlgs
            [~, idxMax] = max(ampEst(:,mc,k0,ii));
            velEst(mc,k0,ii) = velVec(idxMax);
        end
    end
end

%% --- Compute MSE ---
mseVel = mean((velEst - velTrue).^2, 1);
timeTaken_min = (cputime - startTimeTot) / 60;

%% --- Plot Results ---
if doPlot
    figure;
    hold on; grid on;
    title(sprintf('Blood Velocity Estimation (OW = %d), True Velocity = %.2f m/s', no_filt, velTrue));
    xlabel('SNR [dB]'); ylabel('MSE [dB]');

    methodIdxToPlot = [1 2 3 4 5];
    myLegend = {};
    for ii = methodIdxToPlot
        plot(SNR, 10*log10(squeeze(mseVel(1, :, ii))), 'LineWidth', 1.5);
        myLegend{end+1} = algNames{ii};
    end
    legend(myLegend, 'Location', 'northeast');
end

%% --- Save Data ---
if saveData
    save(saveStr, 'ampEst', 'velEst', 'mseVel', 'SNR', 'MC', 'no_filt', 'algNames', 'velTrue');
end
