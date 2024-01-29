clear;
%% Parameters
f_t = 4000; % Frequency of the Tinnitus in Hz
fs = 44100; % Sampling frequency (in Hz)
cycle_dur = 1 / 1.5; % Duration of one cycle in seconds
% cycle_dur = 5; % Duration of one cycle in seconds
protocol_dur = 60; % Protocol duration in s
on2off = [3, 2];
n_side = 2;

% mode = "reg aCR";
mode = "noisy aCR";
% mode = "log";

%% obtain remaining frequencies -> vector
if mode ~= "log"

    S = [0.766, 0.69, 0.728;...
         0.9, 0.81, 0.855;...
         1.1, 1.182, 1.265;...
         1.4, 1.505, 1.61]';
    
    f = S(1, :) * f_t;
else % log mode

    log_range = 0.1; % log distance to mutual frequency components
    exponent1 = log10(f_t) - log_range;
    exponent2 = log10(f_t) + log_range;
    
    vec1 = logspace(exponent1, log10(f_t), n_side+1);
    vec2 = logspace(log10(f_t), exponent2, n_side+1);
    
    f = [vec1(1:end-1), vec2(2:end)];
end

%% Create protocol
n_freqs = length(f);
no_of_cycles = ceil(protocol_dur / cycle_dur);
protocol = zeros(1, n_freqs*no_of_cycles);
for i = 1:no_of_cycles
    if mod(i-1, sum(on2off)) < on2off(1)
        protocol(1+(i-1)*n_freqs : i*n_freqs) = randperm(n_freqs);
    else
        protocol(1+(i-1)*n_freqs : i*n_freqs) = zeros(1, n_freqs);
    end
end

%% Generate the time vector and window function
subcycle_dur = cycle_dur / n_freqs;
t = linspace(0, subcycle_dur, fs * subcycle_dur);
w = hann(subcycle_dur * fs);
window = repmat(w', n_freqs, 1);

%% Generate the signal
% signal = zeros(1, fs * cycle_dur * no_of_cycles);
signal = [];
sin_wave = sin(2 * pi * f' * t) .* window;
zero_wave = zeros(1, subcycle_dur * fs);
for i = 1:length(protocol)
    if(protocol(i) > 0)
%         signal(1+(i-1)*subcycle_dur*fs : i*subcycle_dur*fs) = sin_wave(protocol(i), :);
        signal = [signal, sin_wave(protocol(i), :)];

        if mode == "noisy aCR" % Obtain new f vector
            rand_ids = randi([1, size(S, 1)], 1, size(S, 2));
            rand_lin_ids = rand_ids + ((1:size(S, 2)) - 1) * size(S, 1); % sub2ind
            f_new = S(rand_lin_ids) * f_t;

            sin_wave = sin(2 * pi * f_new' * t) .* window;
        end
    else
        signal = [signal, zero_wave];
    end
end


%% Normalize the signal to be between -1 and 1 and create time vector
signal = signal / max(abs(signal));
t_tot = (1:length(signal)) / fs;

%% Plot the generated signal
figure;
plot(t_tot, signal);
title('Generated time-domain signal');
xlabel('Time (s)');
ylabel('Normalized Amplitude');

%% Save the sine wave as a WAV file
audiowrite('tinnitus_cr.wav', signal, fs);

disp('Sine wave saved as "sine_wave.wav"');
