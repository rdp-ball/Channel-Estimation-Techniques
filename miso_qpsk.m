clc;
clear;

N = 1000; % Number of symbols
x = randi([0 1], 1, 2*N);

% QPSK Modulation (Baseband)
xmod = ((1-2*x(1:2:end)) + 1j * (1-2*x(2:2:end))) / sqrt(2);

% Alamouti STBC Encoding
s1 = xmod(1:2:end);
s2 = xmod(2:2:end);
s = [s1; s2];
s = [s; [-conj(s(2,:)); conj(s(1,:))]];

% Channel: Assume two independent channels (MISO)
h1 = (randn(1, N/2) + 1j * randn(1, N/2)) / sqrt(2);
h2 = (randn(1, N/2) + 1j * randn(1, N/2)) / sqrt(2);

% Noise generation
noise = (randn(2, N/2) + 1j * randn(2, N/2)) / sqrt(2);

snr_db = 0:2:20; % SNR values in dB
ber = zeros(size(snr_db)); % To store BER values

for i = 1:length(snr_db)
    snrdb = snr_db(i);
    snrlin = db2pow(snrdb);

    % Received signal with noise through the channels
    yrx1 = h1 .* s(1,:) + h2 .* s(2,:) + sqrt(1/snrlin) * noise(1,:);
    yrx2 = h1 .* s(3,:) + h2 .* s(4,:) + sqrt(1/snrlin) * noise(2,:);

    % Alamouti STBC Decoding
    y1 = yrx1;
    y2 = yrx2;

    % Combined received signal
    r1 = conj(h1) .* y1 + h2 .* conj(y2);
    r2 = conj(h2) .* y1 - h1 .* conj(y2);

    % Normalize received signals
    h_combined = abs(h1).^2 + abs(h2).^2;
    r1 = r1 ./ h_combined;
    r2 = r2 ./ h_combined;

    % Combine received symbols
    r_combined = reshape([r1; r2], 1, []);

    % QPSK Demodulation
    ydemod = zeros(1, 2*N);
    ydemod(1:4:end) = real(r1) < 0; % In-phase component
    ydemod(2:4:end) = imag(r1) < 0; % Quadrature component
    ydemod(3:4:end) = real(r2) < 0; % In-phase component
    ydemod(4:4:end) = imag(r2) < 0; % Quadrature component

    % BER Calculation
    BERcount = sum(xor(ydemod, x));
    ber(i) = BERcount / (2 * N);
end

% Plot BER vs SNR
figure;
semilogy(snr_db, ber, '-o');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for MISO QPSK with Alamouti STBC');
