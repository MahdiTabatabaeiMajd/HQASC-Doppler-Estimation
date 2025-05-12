function [P_Welch, P_Welch_REF_ham] = REF_Spectral_Ham(Y, N, M)
% -------------------------------------------------------------------------
% REF_Spectral_Ham
%
% Authors:
%   Seyed Mohammad Mahdi Tabatabaei Majd
%   Leila Eslami
%
% Description:
%   This function computes power spectral estimates using Welch's method
%   based on uniform rectangular and Hamming windowing. It calculates the
%   sample covariance matrices of temporally overlapping subarrays with and
%   without Hamming weighting, and estimates the power spectrum over a set
%   of Doppler frequencies.
%
% Inputs:
%   Y : Input data matrix (depth samples x ensemble frames)
%   N : Subarray length (observation window size)
%   M : Number of frequency bins
%
% Outputs:
%   P_Welch         : Power spectrum estimated with rectangular window
%   P_Welch_REF_ham : Power spectrum estimated with Hamming window
%
% -------------------------------------------------------------------------

% Ensure input is oriented as [samples x frames]
Y = transpose(Y);
[Averaging_kernel, K] = size(Y);
L = Averaging_kernel - N + 1;  % Number of subarrays for averaging

% Initialize covariance matrices
R = zeros(N, N);       % For rectangular window
R_ham = zeros(N, N);   % For Hamming window

% Loop over L overlapping segments to compute average covariance
for l = 1:L
    Y_seg = Y(l:l+N-1, :);  % Extract subarray

    % Rectangular window (no weighting)
    R = R + Y_seg * Y_seg';

    % Hamming window applied across depth (normalized)
    w_ham = repmat(hamming(N), 1, size(Y_seg, 2)) / sum(hamming(N));
    Y_seg_ham = Y_seg .* w_ham;
    R_ham = R_ham + Y_seg_ham * Y_seg_ham';
end

% Normalize covariance matrices
R = R / (K * L);
R_ham = R_ham / (K * L);

% Frequency vector for Doppler FFT
omega_si = 2 * pi * ([0:M-1]' / M - 0.5);
n = (0:N-1)';

% Initialize outputs
P_Welch = zeros(M, 1);
P_Welch_REF_ham = zeros(M, 1);

% Calculate power spectrum at each frequency bin
for no = 1:M
    e = exp(1i * omega_si(no) * n);  % Steering vector

    % Welch with rectangular window
    P_Welch(no) = abs(e' * R * e);

    % Welch with Hamming window
    P_Welch_REF_ham(no) = abs(e' * R_ham * e);
end

end
