% Parameters
N = 64;                % Number of subcarriers
numSymbols = 100;      % Number of OFDM symbols
L = 16;                % Length of the channel (number of taps)
SNR_dB = 20;           % Signal to Noise Ratio in dB
pilotSpacing = 4;      % Pilot spacing

% Generate random bits and map to QPSK symbols
bitsPerSymbol = 2;     % QPSK has 2 bits per symbol
totalBits = numSymbols * N * bitsPerSymbol;
bits = randi([0 1], totalBits, 1);

% QPSK Modulation
symbols = 1/sqrt(2) * ((1-2*bits(1:2:end)) + 1i*(1-2*bits(2:2:end)));

% Reshape symbols into OFDM symbols
txSymbols = reshape(symbols, N, numSymbols);

% Print Input Symbols
fprintf('Input Symbols (First OFDM symbol):\n');
disp(txSymbols(:, 1));

% IFFT (OFDM modulation)
txOFDM = ifft(txSymbols, N);

% Add cyclic prefix
cyclicPrefixLength = L - 1;
txOFDM_CP = [txOFDM(end-cyclicPrefixLength+1:end, :); txOFDM];

% Create channel (Rayleigh fading channel)
h = (1/sqrt(2)) * (randn(L, 1) + 1i*randn(L, 1));

% Print the Channel
fprintf('Channel (Impulse Response):\n');
disp(h);

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

% LS Channel Estimation
% Extract pilot symbols
pilotIndices = 1:pilotSpacing:N;
txPilots = txSymbols(pilotIndices, :);
rxPilots = rxSymbols(pilotIndices, :);

% Print Pilot Symbols
fprintf('Pilot Symbols (Transmitted):\n');
disp(txPilots(:, 1));

fprintf('Pilot Symbols (Received):\n');
disp(rxPilots(:, 1));

% LS estimation of channel on pilot positions
H_ls_pilot = rxPilots ./ txPilots;

% Print LS estimation at Pilot positions
fprintf('LS Estimation at Pilot Positions (First OFDM symbol):\n');
disp(H_ls_pilot(:, 1));

% Interpolate the channel estimate across all subcarriers
H_ls = interp1(pilotIndices, H_ls_pilot, 1:N, 'linear', 'extrap');

% Print Final Channel Estimation
fprintf('Final Channel Estimation (LS) for all Subcarriers (First OFDM symbol):\n');
disp(H_ls.');

% Equalization
equalizedSymbols = rxSymbols ./ H_ls;

% Demapping
receivedBits = zeros(totalBits, 1);
receivedBits(1:2:end) = real(equalizedSymbols(:)) < 0;
receivedBits(2:2:end) = imag(equalizedSymbols(:)) < 0;

% Calculate Bit Error Rate (BER)
BER = sum(bits ~= receivedBits) / length(bits);

% Display the results
fprintf('SNR: %d dB\n', SNR_dB);
fprintf('BER: %e\n', BER);

% Plot estimated and actual channel response
figure;
subplot(3,1,1);
plot(1:N, abs(fft(h, N)), 'b', 'LineWidth', 1.5);
title('Actual Channel Frequency Response');
xlabel('Subcarrier Index');
ylabel('Magnitude');
grid on;

subplot(3,1,2);
plot(1:N, abs(H_ls), 'r', 'LineWidth', 1.5);
title('Estimated Channel Frequency Response (LS)');
xlabel('Subcarrier Index');
ylabel('Magnitude');
grid on;

% Plot Pilot Symbols
subplot(3,1,3);
stem(pilotIndices, abs(txPilots(:, 1)), 'bo', 'filled', 'DisplayName', 'Transmitted Pilots');
hold on;
stem(pilotIndices, abs(rxPilots(:, 1)), 'rx', 'DisplayName', 'Received Pilots');
hold off;
title('Pilot Symbols');
xlabel('Pilot Subcarrier Index');
ylabel('Magnitude');
legend('show');
grid on;
