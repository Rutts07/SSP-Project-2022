close all;
clc;

% Load the speech signal
% data -> number of samples
% fs -> sampling frequency
[data, fs] = audioread('sa2.wav');
% sound(data, fs);

% Apply framing onto the signal
frame_size = 0.02;                                  % 20ms frames
frame_over = 0.5;                                   % 50% overlap
frame_samp = fs*frame_size;                         % samples in the frame
num_frames = floor(length(data)/frame_samp)*2 - 1;  % number of frames
spec_peaks = zeros(num_frames, 1);

i = 1; j = 2;
while i < length(data)
    lower = i;
    upper = i + frame_samp;
    
    if upper > length(data)
        upper = length(data);
    end

    frame_sign = data(lower:upper);
    % frame_sign = frame_sign*hamming(length(frame_sign));

    frame_nfft = fft(frame_sign, 256);
    frame_nfft = frame_nfft(1:128);
    [sorted, index] = sort(frame_nfft, 'descend');
    max10 = sorted(1:10);
    
    amp = 0;
    for temp=1:10
        amp = amp + abs(max10(temp));
    end

    spec_peaks(j) = amp;

    i = i + (frame_samp*frame_over);
    j = j + 1;
end

spec_peaks = smooth(spec_peaks);

% Normalize the FOD spectral peaks
norm_peaks = spec_peaks;
enhanced_spec_peaks = zeros(num_frames, 1);
norm_samps = 0;
slope = 1;
lower = 20;
upper = 20;

while norm_samps < num_frames
    if norm_peaks(upper) < norm_peaks(upper + 1)
        slope = 1;
    elseif norm_peaks(upper) > norm_peaks(upper + 1)
        slope = 0;
    else
        slope = -1;
        upper = upper + 1;
    end

    % positive slope
    if slope == 1 && upper < num_frames
        while norm_peaks(upper) < norm_peaks(upper + 1)
            upper = upper + 1;
            if upper >= num_frames
                break;
            end
        end

        if norm_peaks(upper) - norm_peaks(lower) > 10
            for j=lower:upper
                if j <= num_frames
                    enhanced_spec_peaks(j, 1) = (norm_peaks(j) - norm_peaks(lower)) / (norm_peaks(upper) - norm_peaks(lower));            
                end
            end
        else
            for j=lower:upper
                if j <= num_frames
                    enhanced_spec_peaks(j, 1) = enhanced_spec_peaks(lower-1);
                end
            end
        end

        lower = upper + 1;
        norm_samps = upper + 1;
    end

    if norm_peaks(upper) < norm_peaks(upper + 1)
        slope = 1;
    elseif norm_peaks(upper) > norm_peaks(upper + 1)
        slope = 0;
    else
        slope = -1;
        upper = upper + 1;
    end

    % negative slope
    if slope == 0 && upper < num_frames
        while norm_peaks(upper) > norm_peaks(upper + 1)
            upper = upper + 1;
            if upper >= num_frames
                break;
            end
        end

        if norm_peaks(lower) - norm_peaks(upper) > 10
            for j=lower:upper
                if j <= num_frames
                    enhanced_spec_peaks(j, 1) = (norm_peaks(j) - norm_peaks(upper)) / (norm_peaks(lower) - norm_peaks(upper));
                end
            end
        else
            for j=lower:upper
                if j <= num_frames
                    enhanced_spec_peaks(j, 1) = enhanced_spec_peaks(lower-1);
                end
            end
        end

        lower = upper + 1;
        norm_samps = upper + 1;
    end
end
enhanced_spec_peaks = smooth(enhanced_spec_peaks);

%{
% Enhance the spectral peaks with smoothed hilbert transform and then FOD
enhanced_spec_peaks = smooth(diff(smooth(hilbert(enhanced_spec_peaks))));
%}

% Enhanced spectrum sent to FOGD
fogd_size = 0.1*fs;                     % 100ms windows
wind_size = fogd_size / frame_samp;     % num_frames to convolve
gaussFilter = gausswin(wind_size);
fogd_spectrum = conv(gaussFilter, enhanced_spec_peaks);
fogd_spectrum = smooth(fogd_spectrum);

% Apply smoothing several times x20
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);
fogd_spectrum = smooth(fogd_spectrum);

% Create a time axis
time_sign = linspace(0, length(data)/fs, length(data));
time_fram = linspace(0, num_frames*frame_size*frame_over, length(spec_peaks));
time_fogd = linspace(0, num_frames*frame_size*frame_over, length(fogd_spectrum));

figure;
subplot(5, 1, 1);
plot(time_sign, data);
title('Time plot of speech signal');

subplot(5, 1, 2);
plot(linspace(0, num_frames*frame_size*frame_over, length(spec_peaks)), spec_peaks);
title('Sum of ten largest peaks in the DFT spectrum');

subplot(5, 1, 3);
plot(linspace(0, num_frames*frame_size*frame_over, length(enhanced_spec_peaks)), real(enhanced_spec_peaks));
title('Enhanced sum of ten largest peaks in the DFT spectrum');

subplot(5, 1, 4);
x = time_fogd; 
y = normalize(real(fogd_spectrum));
y = y - min(y);
y = y/max(y);

% Vowel Onset Points
% pos_peaks = islocalmax(y);
[pos_peaks, pos_peaks_idx] = findpeaks(y);

%{
% Eliminate Spurious Peaks
pos_peaks = pos_peaks(3:end-2);
pos_peaks_idx = pos_peaks_idx(3:end-2);
%}

% Vowel Offset Points
% neg_peaks = islocalmin(y);

plot(x, y);
hold on
stem(x(pos_peaks_idx), pos_peaks, '*');
% plot(pos_peaks, y(pos_peaks_idx), '^r', 'Color', [1 0 0], 'MarkerSize', 5);
% plot(x(neg_peaks), y(neg_peaks), 'vr', 'Color', [0 0 1], 'MarkerSize', 2);
hold off
title('VOP Evidence Plot');

subplot(5, 1, 5);
plot(time_sign, data);
title('Time plot of speech signal');
hold on
stem(x(pos_peaks_idx), pos_peaks/15, '*');
hold off
title('Vowel Onset Points');
