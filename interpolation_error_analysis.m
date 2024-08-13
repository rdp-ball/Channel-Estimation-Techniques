clc;clear;close all;
interpolation_error_analysis('linear');
interpolation_error_analysis('spline');
interpolation_error_analysis('cubic');
interpolation_error_analysis('pchip');

function interpolation_error_analysis(interpolation_method)
    t = 0:0.1:2*pi; 
    %original_signal = sin(5*t) + 0.5*sin(10*t);
    
    original_signal=sin(2*pi*5*t) + 0.5*sin(2*pi*10*t) + 0.3*sin(2*pi*20*t) + 0.2*sin(2*pi*30*t);
    
    % Sample the original signal after 2 points
    sampled_points = original_signal(1:2:end); 

    % Do Interpolation
    t_interp = 0:0.1:2*pi; 
    interp_signal = interp1(t(1:2:end), sampled_points, t_interp, interpolation_method); 

    % FInding Error 
    difference = original_signal - interp_signal;
    nan_indices = isnan(interp_signal);
    difference = difference(~nan_indices);

   
    if norm(original_signal) ~= 0 && norm(interp_signal) ~= 0
        mean_abs_error = mean(abs(difference));
        normalized_error = mean_abs_error / norm(original_signal);
        disp(['Normalized Error between original and interpolated signals (', interpolation_method, '): ', num2str(normalized_error)]);
        disp(['Mean Absolute Error between original and interpolated signals (', interpolation_method, '): ', num2str(mean_abs_error)]);
    else
        disp(['Unable to calculate mean absolute error. Original or interpolated signal has zero norm. (', interpolation_method, ')']);
        return;
    end

    % Plotting the signals
    figure;
    subplot(3,1,1);
    plot(t, original_signal, 'b-', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Amplitude');
    title('Original Signal');
    grid on;

    subplot(3,1,2);
    stem(t(1:2:end), sampled_points, 'k', 'filled');
    xlabel('Time');
    ylabel('Amplitude');
    title('Sampled Points');
    grid on;

    subplot(3,1,3);
    plot(t_interp, interp_signal, 'r--', 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Amplitude');
    title('Interpolated Signal');
    grid on;

    % Comparison of Original Signal vs Interpolated Signal
    figure;
    plot(t, original_signal, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_interp, interp_signal, 'r--', 'LineWidth', 1.5);
    scatter(t(1:2:end), sampled_points, 'k', 'filled');
    legend('Original Signal', 'Interpolated Signal', 'Sampled Points');
    xlabel('Time');
    ylabel('Amplitude');
    title(['Comparison of Original and Interpolated Signals (', interpolation_method, ')']);
    grid on;
    hold off;
end
