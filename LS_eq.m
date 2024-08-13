clc
N = 10;
% Generate random BPSK symbols
symbols = randi([0 1], 1, N)
% Map BPSK symbols to complex constellation points
mapping_table = exp(1j * pi * [0 1]); % BPSK constellation: 0 maps to 1, 1 maps to -1
x_mod = mapping_table(symbols + 1); % Mapping symbols to BPSK constellation points
h = (randn(1, N) + 1j * randn(1, N)) / sqrt(2); % Channel impulse response (Rayleigh fading channel)
n = randn(1, N) + 1j * randn(1, N); % Additive White Gaussian Noise (AWGN) with zero mean and unit variance
yrx = h .* x_mod + 0.1 * n; % Received signal after passing through the channel and adding noise
scatterplot(yrx./h); % Scatter plot of received signal after equalization

% Calculate SNR
signal_power = mean(abs(h .* x_mod).^2); % Power of the transmitted signal
noise_power = mean(abs(n).^2); % Power of the noise
SNR = signal_power / noise_power; % Signal-to-Noise Ratio
disp(['SNR: ', num2str(10*log10(SNR)), ' dB']);

% LS Equalizer
y_demod = yrx ./ h; % Equalization by division
% Demodulation (BPSK)
demod_symbols = real(y_demod) < 0; % Demodulate by comparing the real part of y_demod with 0
BER_count = sum(demod_symbols ~= symbols); % Count Bit Errors for BPSK
disp(['Bit Error Rate (BER): ', num2str(BER_count / N)]);

% Display the demodulated signal
disp('Demodulated Symbols:');
disp(demod_symbols);
