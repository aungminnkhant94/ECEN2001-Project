% Given parameters
N = 128; % Length of ZC sequence
N_CP = 16; % Length of cyclic prefix
q = 3; % ZC sequence parameter

% Generate Zadoff–Chu sequence
n = (0:N-1)';
ZC_sequence = exp(-1i * pi * q * n .* (n + mod(N, 2)) / N);

% Add cyclic prefix to the ZC sequence
ZC_with_prefix = [ZC_sequence(end - N_CP + 1:end); ZC_sequence];

% Example: Generating a random received signal for demonstration
received_signal_length = 1000; % Length of received signal (adjust as needed)
y = randn(1, received_signal_length) + 1i * randn(1, received_signal_length); % Example of received signal with AWGN

% Perform cross-correlation to detect preamble
correlation_result = abs(xcorr(1,256));

% Find the peaks in the cross-correlation result
[peaks, peak_indices] = findpeaks(correlation_result, 'SortStr', 'ascend');

% Find the highest peak
if ~isempty(peaks)
max_peak = peaks(1);
max_peak_index = peak_indices(1);
disp(['Highest Peak Value: ', num2str(max_peak)]);
disp(['Index of Highest Peak: ', num2str(max_peak_index)]);
else
disp('No peaks found.');
end

% Find the second peak if available
if length(peaks) > 1
    % Exclude the index of the highest peak found earlier
    peak_indices_without_max = peak_indices(peak_indices ~= max_peak_index);
    
    % Find the second highest peak among the remaining peaks
    second_max_peak = peaks(peak_indices_without_max == max(peak_indices_without_max));
    second_max_peak_index = peak_indices_without_max(peak_indices_without_max == max(peak_indices_without_max));
    
    disp(['Second Highest Peak Value: ', num2str(second_max_peak)]);
    disp(['Index of Second Highest Peak: ', num2str(second_max_peak_index)]);
end

% Plot the cross-correlation result with detected peaks
figure;
plot(correlation_result);
hold on;
plot(peak_indices, peaks, 'r*', 'MarkerSize', 10);
xlabel('Sample Index');
ylabel('Correlation Value');
title('Cross-Correlation with ZC Sequence');
legend('Cross-Correlation Result', 'Detected Peaks');
hold off;

% Get the index of the maximum peak
[max_peak, max_peak_index] = max(peaks);
disp(['Maximum Peak Value: ', num2str(max_peak)]);
disp(['Index of Maximum Peak: ', num2str(peak_indices(max_peak_index))]);


% Peak Detection
if ~isempty(peaks)
    % Determine if the highest peak meets a certain threshold
    peak_threshold = 0.5 * max_peak; % Adjust threshold as needed
    if max_peak >= peak_threshold
        disp('ZC sequence detected.'); % ZC sequence detected
    else
        disp('No ZC sequence detected.'); % ZC sequence not detected
    end
else
    disp('No peaks found.'); % No peaks found
end

% Timing Synchronization
if ~isempty(peaks)
    % Determine timing offset based on the index of the highest peak
    timing_offset = max_peak_index - N_CP; % Timing offset in samples
    disp(['Timing Offset: ', num2str(timing_offset)]);
else
    disp('No peaks found.'); % No peaks found
end

% Frame Detection
if ~isempty(peaks) && (max_peak >= peak_threshold)
    % Remove cyclic prefix and extract data portion of received signal
    data_without_prefix = y(max_peak_index:end); % Data portion after preamble
    % Further processing of the data...
else
    disp('Cannot proceed with frame detection.'); % Peak not above threshold or no peaks found
end

% Mapping steps
% Given parameters
N = 128; % Number of subcarriers
M = 16; % QAM modulation order
num_symbols = N; % Number of QAM symbols, same as the number of subcarriers

% Generate random QAM symbols (alternative approach)
qam_symbols = qammod(randi([0 M-1], num_symbols, 1), M); % Random QAM symbols

% Perform Discrete Fourier Transform (DFT) to map QAM symbols to frequency domain
frequency_domain_symbols = fft(qam_symbols);

% Convert to time-domain waveform using Inverse Discrete Fourier Transform (IDFT)
time_domain_waveform = ifft(frequency_domain_symbols);

% Normalize the time-domain waveform (optional)
time_domain_waveform = time_domain_waveform / max(abs(time_domain_waveform));

% Plot the time-domain waveform (optional)
figure;
plot(abs(time_domain_waveform));
xlabel('Sample Index');
ylabel('Amplitude');
title('Time-Domain Waveform');
