function MIMO_Bayesian_Channel_Estimation()
    % System parameters
    Nt = 3; % Number of transmit antennas
    Nr = 4; % Number of receive antennas
    Np = 10; % Number of pilot symbols
    Nd = 100; % Number of data symbols
    SNR_dB = 0:5:30; % SNR range in dB
    num_iterations = 1000; % Number of iterations 

    % Pre-allocate BER array
    BER = zeros(length(SNR_dB), 1);

    % Loop over different SNR values
    for idx = 1:length(SNR_dB)
        SNR = 10^(SNR_dB(idx) / 10);
        sigma2_n = 1 / SNR;

        
        error_count = 0;
        total_bits = 0;

        for iteration = 1:num_iterations
            % Generate random pilot signals
            X = sqrt(1) * (randn(Nt, Np) + 1i * randn(Nt, Np));

            % Generate a random channel matrix
            H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);

            Y_pilot = H * X + sqrt(sigma2_n / 2) * (randn(Nr, Np) + 1i * randn(Nr, Np));

            
            H_est = BayesianChannelEstimation(Y_pilot, X, sigma2_n, Nt, Nr);

            % Generate random data symbols
            data_symbols = randi([0 1], Nt, Nd);
            data_symbols_mod = 2 * data_symbols - 1; % BPSK modulation

            % Combine pilot and data symbols
            X_combined = [X, data_symbols_mod];

            % Transmit combined symbols through the channel
            Y_combined = H * X_combined + sqrt(sigma2_n / 2) * (randn(Nr, Np+Nd) + 1i * randn(Nr, Np+Nd));

            % Separate received pilot and data signals
            Y_data = Y_combined(:, Np+1:end);

            % Receiver side: 
            data_symbols_est = pinv(H_est) * Y_data;
            data_symbols_est_demod = real(data_symbols_est) > 0;

            % Count bit errors
            error_count = error_count + sum(sum(data_symbols ~= data_symbols_est_demod));
            total_bits = total_bits + numel(data_symbols);
        end

        % Compute BER
        BER(idx) = error_count / total_bits;
    end

    % Plot BER vs SNR
    figure;
    semilogy(SNR_dB, BER, 'o-');
    xlabel('SNR (dB)');
    ylabel('Bit Error Rate (BER)');
    title('BER vs SNR for MIMO Bayesian Channel Estimation');
    grid on;

    disp('True Channel Matrix:');
    disp(H);

    disp('Estimated Channel Matrix:');
    disp(H_est);

    disp('Original Data Symbols:');
    disp(data_symbols);

    disp('Pilot Signals:');
    disp(X);
end

function H_est = BayesianChannelEstimation(Y, X, sigma2_n, Nt, Nr)
    P = 2; % Transmit power

    % Compute the posterior mean and covariance
    Sigma_H = inv((1 / sigma2_n) * (X * X') + (1 / P) * eye(Nt));
    H_mean = (1 / sigma2_n) * Sigma_H * X * Y';

    
    H_est = H_mean';

   
end
