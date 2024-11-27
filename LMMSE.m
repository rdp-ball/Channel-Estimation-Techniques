clc;clear;close all;
% Parameters
N = 64;                % Number of subcarriers
L = 16;                % Length of the channel (number of taps)
pilotSpacing = 4;      % Pilot spacing
SNR_range = 0:5:30;    % Range of SNR values in dB for BER vs. SNR plot
numSNR = length(SNR_range);  % Number of SNR points
numSymbols = 1000;     % Number of OFDM symbols
numIterations = 50;    % Number of iterations to average results

% Initialize BER array
BER_LMMSE = zeros(1, numSNR);

% Loop over different SNR values
for idx = 1:numSNR
    SNR_dB = SNR_range(idx);
    berSum = 0;  % Accumulator for BER over iterations
    
    for iter = 1:numIterations
        % Generate random bits and map to QPSK symbols
        bitsPerSymbol = 2;     % QPSK has 2 bits per symbol
        totalBits = numSymbols * N * bitsPerSymbol;
        bits = randi([0 1], totalBits, 1);

        % QPSK Modulation
        symbols = 1/sqrt(2) * ((1-2*bits(1:2:end)) + 1i*(1-2*bits(2:2:end)));

        % Reshape symbols into OFDM symbols
        txSymbols = reshape(symbols, N, numSymbols);

        % IFFT (OFDM modulation)
        txOFDM = ifft(txSymbols, N);

        % Add cyclic prefix
        cyclicPrefixLength = L - 1;
        txOFDM_CP = [txOFDM(end-cyclicPrefixLength+1:end, :); txOFDM];

        % Create channel (Rayleigh fading channel)
        h = (1/sqrt(2)) * (randn(L, 1) + 1i*randn(L, 1));

        % Convolve transmitted signal with channel
        rxOFDM_CP = filter(h, 1, txOFDM_CP);

        % Add noise
        noisePower = 10^(-SNR_dB/10);
        noise = sqrt(noisePower/2) * (randn(size(rxOFDM_CP)) + 1i*randn(size(rxOFDM_CP)));
        rxOFDM_CP = rxOFDM_CP + noise;

        % Remove cyclic prefix
        rxOFDM = rxOFDM_CP(cyclicPrefixLength+1:end, :);

        % FFT (OFDM demodulation)
        rxSymbols = fft(rxOFDM, N);

        % LS Channel Estimation (used as initial estimate)
        pilotIndices = 1:pilotSpacing:N;
        txPilots = txSymbols(pilotIndices, :);
        rxPilots = rxSymbols(pilotIndices, :);

        H_ls_pilot = rxPilots ./ txPilots;

        % LMMSE Channel Estimation
        % Noise Variance
        sigma_n2 = noisePower;

        % Channel correlation matrix (assuming exponential decay model)
        rho = 0.9; % Correlation coefficient
        Rhh = toeplitz(rho.^(0:N-1));

        % Estimate noise power spectral density
        SNR_linear = 10^(SNR_dB/10);
        sigma_h2 = var(h);  % Channel power
        sigma_n2 = sigma_h2 / SNR_linear;

        % LMMSE filter
        H_lmmse = zeros(N, numSymbols);

        for i = 1:numSymbols
            % Extract LS estimates for current symbol
            H_ls = interp1(pilotIndices, H_ls_pilot(:, i), 1:N, 'linear', 'extrap').'; % Ensure column vector

            % Calculate LMMSE estimate
            R_inv = Rhh * inv(Rhh + (sigma_n2/sigma_h2)*eye(N));  % Updated LMMSE filter expression
            H_lmmse(:, i) = R_inv * H_ls;
        end

        % Equalization
        equalizedSymbols = rxSymbols ./ H_lmmse;

        % Demapping
        receivedBits = zeros(totalBits, 1);
        receivedBits(1:2:end) = real(equalizedSymbols(:)) < 0;
        receivedBits(2:2:end) = imag(equalizedSymbols(:)) < 0;

        % Calculate Bit Error Rate (BER)
        berSum = berSum + sum(bits ~= receivedBits) / length(bits);
    end
    
    % Average BER over iterations
    BER_LMMSE(idx) = berSum / numIterations;
end

% Plot BER vs. SNR
figure;
semilogy(SNR_range, BER_LMMSE, 'bo-', 'LineWidth', 1.5);
title('BER vs. SNR (LMMSE)');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
grid on;

% Plot estimated and actual channel response for one OFDM symbol
figure;
subplot(2,1,1);
plot(1:N, abs(fft(h, N)), 'b', 'LineWidth', 1.5);
title('Actual Channel Frequency Response');
xlabel('Subcarrier Index');
ylabel('Magnitude');
grid on;

subplot(2,1,2);
plot(1:N, abs(H_lmmse(:, 1)), 'r', 'LineWidth', 1.5);
title('Estimated Channel Frequency Response (LMMSE)');
xlabel('Subcarrier Index');
ylabel('Magnitude');
grid on;

% Display the results
fprintf('SNR: %d dB\n', SNR_dB);
fprintf('BER: %e\n', BER_LMMSE(end));

% Print the first OFDM symbol's input and output symbols and the channel
fprintf('Input Symbols (First OFDM symbol):\n');
disp(txSymbols(:, 1));

fprintf('Output Symbols (First OFDM symbol):\n');
disp(rxSymbols(:, 1));

fprintf('Channel (Impulse Response):\n');
disp(h);
