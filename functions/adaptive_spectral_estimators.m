function [P_Welch, P_Capon, P_P_C, MASC, HQASC] = adaptive_spectral_estimators(Y, N, M, W_t)
% -------------------------------------------------------------------------
% adaptive_spectral_estimators
%
% Authors:
%   Seyed Mohammad Mahdi Tabatabaei Majd
%   Leila Eslami
%
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
%
% Description:
%   This function estimates Doppler power spectra using several spectral
%   estimation techniques, including Welch, Capon, Projected Capon (Pr.Capon),
%   MASC, and HQASC. The method uses sample covariance estimation
%   and adaptive eigen-decomposition for high-resolution Doppler imaging.
%
% Inputs:
%   Y   : Input data matrix (samples x ensembles)
%   N   : Sub-aperture size
%   M   : FFT size (number of frequency bins)
%   W_t : Temporal frequency index
%
% Outputs:
%   P_Welch  : Power spectrum estimated by Welch method
%   P_Capon  : Power spectrum estimated by Capon method
%   P_P_C    : Projected Capon power spectrum
%   MASC     : Eigen-based adaptive spectral coherence estimator
%   HQASC    : High-resolution adaptive spectral coherence (block 1)
% -------------------------------------------------------------------------

Y = transpose(Y);  % Ensure data is in [samples x ensembles] form
[Averaging_kernel, K] = size(Y);
L = Averaging_kernel - N + 1;  % Number of sliding windows
J = fliplr(eye(N));            % Anti-diagonal identity matrix for FB averaging

% --- Estimate sample covariance matrix R ---
R = zeros(N, N);
for l = 1:L
    Y_zir = Y(l:l+N-1, :);
    R = R + Y_zir * Y_zir';
end
R = R / (K * L);

% Forward-backward averaged covariance matrix
RFB = (R + J * R.' * J) / 2;

% --- Eigen decomposition and separation ---
[eigvectors, eigvalus] = eig(R);
DD = fliplr(flipud(eigvalus));
lamda_max = max(DD(:));
num = 0;
for p = 1:size(DD,1)
    if DD(p,p) >= 0.5 * lamda_max
        num(p) = DD(p,p);
    end
end

V = fliplr(eigvectors);
Es = V(:,1:nnz(num));      % Signal subspace
En = V(:,nnz(num)+1:end);  % Noise subspace
R_inv = inv(R);

% --- FB averaged eigen decomposition ---
[eigvectorsFB, eigvalusFB] = eig(RFB);
DDFB = fliplr(flipud(eigvalusFB));
lamda_maxFB = max(DDFB(:));
numFB = 0;
for pFB = 1:size(DDFB,1)
    if DDFB(pFB,pFB) >= 0.5 * lamda_maxFB
        numFB(pFB) = DDFB(pFB,pFB);
    end
end

VFB = fliplr(eigvectorsFB);
EsFB = VFB(:,1:nnz(numFB));
EnFB = VFB(:,nnz(numFB)+1:end);
R_invFB = inv(RFB);

% Frequency vector
omega_si = 2 * pi * ([0:M-1]' / M - 0.5);
m = [0:N-1]';

% Initialize output spectra
P_Welch = zeros(M,1);
P_Capon = zeros(M,1);
P_P_C = zeros(M,1);
MASC = zeros(M,1);
HQASC = zeros(M,1);

for no = 1:length(omega_si)
    e = exp(1i * omega_si(no) * m);

    % ASC and FB ASC weights
    W_ASC = (R_inv * e) / (e' * R_inv * e);
    Ws = Es * Es' * W_ASC;
    W_ASCFB = (R_invFB * e) / (e' * R_invFB * e);
    WsFB = EsFB * EsFB' * W_ASCFB;

    % Coherence signal power
    CS_v1 = W_ASC' * R * W_ASC;

    % Initialize for averaging
    R_stotal = zeros(N,N);

    for k=1:K
        % Forward data
        G{k} = zeros(N,L);
        T{k} = zeros(1,L);
        G_S{k} = zeros(N,L);
        G_SFB{k} = zeros(N,L);
        for q=1:L
            % Normal Output
            G{k}(:,q) = Y(q:q+N-1,k).*exp(-1i*(W_t*(k-1)+omega_si(no)*(q-1)));

            % Modified Output
            g_in11FB{k}(1,q) = WsFB'*(Y(q:q+N-1,k).*exp(-1i*(W_t*(k-1)+omega_si(no)*(q-1))));

            % Total Energy as Subarray
            T{k}(1,q) = (1/N)*(Y(q:q+N-1,k)'*Y(q:q+N-1,k));
        end
        g{k} = (mean(G{k},2));
        IC(k) = mean(T{k},2);

        %% The High Resolution CS for CF
        g_out11FB(k) = mean(g_in11FB{k},2);

        %% Coherence Factors
        % Zhao based CF
        CF_v_MASC(k) = (CS_v1/(CS_v1+(1/N)*(IC(k)-CS_v1)));         %% Zhao-CF

    end

    %% Modified Nonliear CF With Zhao CF and LI and Digonal Loading

    for kk=1:K

        %% Modified High Resolution Coherence Factor
        % Modifying Block 1
        g_out1FB(kk) = mean(g_in11FB{kk},2);
        g_out2FB(kk) = sqrt(abs(W_ASCFB'*RFB*W_ASCFB));
        out1FB = g_out1FB(kk) + g_out2FB(kk);
        out2FB = g_out1FB(kk) - g_out2FB(kk);
        Out_Modified(kk) = 0.25*abs((abs(out1FB).^2 - abs(out2FB).^2));

        % Modified High Resolution Coherence Factor, HQASC
        CF_vFB_B1(kk) = ((Out_Modified(kk)))/(IC(kk));

    end
    %%  Calculation Power Estimation

    %   Periodogram Estimation
    P_Welch(no) = abs(e'*R*e);

    %   Capon Estimation
    P_Capon(no) = 1/abs(e'*R_inv*e);

    %   Proj.Capon Estimation
    P_P_C(no) =abs(Ws'*R*Ws);



    %%
    MASC(no)=0;
    HQASC(no)=0;

    for ii=1:K

        %   MASC
        MASC(no) = MASC(no)+abs(CF_v_MASC(ii)*Ws'*g{ii}).^2;

        %   HQASC
        HQASC(no) = HQASC(no)+abs(CF_vFB_B1(ii)*WsFB'*g{ii}).^2;
      
    end


end

MASC=MASC/K;
HQASC=HQASC/K;
end
