clc
% Define parameters
N = 100; % Length of training sequence
M = 10; % Length of channel response
SNR_dB = 20; % Signal-to-Noise Ratio (dB)

% Generate random training symbols (known to both transmitter and receiver)
training_symbols = randi([0 1], 1, N);

% Modulate training symbols
mapping_table = exp(1j * pi * [0 1]); % BPSK constellation: 0 maps to 1, 1 maps to -1
x_mod_train = mapping_table(training_symbols + 1); % Modulate training symbols

% Generate random channel response (known to transmitter but unknown to receiver)
h_true = (randn(1, M) + 1j * randn(1, M)) / sqrt(2); % Channel impulse response (Rayleigh fading channel)

% Generate random noise
noise_power = 10^(-SNR_dB/10); % Noise power
n = sqrt(noise_power/2) * (randn(1, N + M - 1) + 1j * randn(1, N + M - 1)); % AWGN

% Transmit signal through channel
y_train = conv(h_true, x_mod_train) + n; % Received signal after passing through the channel and adding noise

% Formulate equation system
X_train = toeplitz([zeros(1, M-1) x_mod_train zeros(1, N-1)], zeros(1, M)); % Data matrix
Y_train = y_train(M:end); % Output vector

% Solve for channel response using LS method
h_est = pinv(X_train) * Y_train'; % Least Squares solution

% Display estimated and true channel response
disp('True Channel Response:');
disp(h_true);
disp('Estimated Channel Response:');
disp(h_est');
