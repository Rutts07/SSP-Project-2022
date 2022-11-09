close all;
clc;

% Create a new directory to store vowel files
train_dir = dir('traindata/*');
mkdir('vowels_train');

for i=3:length(train_dir)   
    % Create directories for each language to store the vowel data
    % train_dir(i) -> language name
    allwav = dir(fullfile('traindata', train_dir(i).name, '*.wav'));
    mkdir(strcat('vowels_train/', train_dir(i).name));
    
    % Loop over all .wav files present in the directory
    for j=1:length(allwav)
        % Save the new vowel file under this name
        org_file = fullfile('traindata', train_dir(i).name, allwav(j).name);
        new_file = fullfile('vowels_train', train_dir(i).name, strcat('vowel_', allwav(j).name));

        % VOP detection code to loop over

        % Load the speech signal
        [data, fs] = audioread(org_file);
        % data -> number of samples, fs -> sampling frequency

        % Apply framing onto the signal
        frame_size = 0.02;                                  % 20ms frames
        frame_over = 0.5;                                   % 50% overlap
        frame_samp = fs*frame_size;                         % samples in the frame
        num_frames = floor(length(data)/frame_samp)*2 - 1;  % number of frames
        spec_peaks = zeros(num_frames, 1);                 

        k = 1; l = 2;
        while k < length(data)
            lower = k;
            upper = k + frame_samp;
    
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

            spec_peaks(l) = amp;

            k = k + (frame_samp * frame_over);
            l = l + 1;
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
                    for l=lower:upper
                        if l <= num_frames
                            enhanced_spec_peaks(l, 1) = (norm_peaks(l) - norm_peaks(lower)) / (norm_peaks(upper) - norm_peaks(lower));            
                        end
                    end
                else
                    for l=lower:upper
                        if l <= num_frames
                            enhanced_spec_peaks(l, 1) = enhanced_spec_peaks(lower-1);
                        end
                    end
                end

                lower = upper + 1;
                norm_samps = upper + 1;
            end

            if upper >= num_frames
                break;
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
                    for l=lower:upper
                        if l <= num_frames
                            enhanced_spec_peaks(l, 1) = (norm_peaks(l) - norm_peaks(upper)) / (norm_peaks(lower) - norm_peaks(upper));
                        end
                    end
                else
                    for l=lower:upper
                        if l <= num_frames
                            enhanced_spec_peaks(l, 1) = enhanced_spec_peaks(lower-1);
                        end
                    end
                end

                lower = upper + 1;
                norm_samps = upper + 1;

            if upper >= num_frames
                break;
            end

            end
        end

        enhanced_spec_peaks = smooth(enhanced_spec_peaks);

        % Enhance the spectral peaks with smoothed hilbert transform
        % enhanced_spec_peaks = smooth(hilbert(enhanced_spec_peaks));

        % Enhanced spectrum sent to FOGD
        fogd_size = 0.1*fs;                     % 100ms windows
        wind_size = fogd_size / frame_samp;     % num_frames to convolve
        gaussFilter = gausswin(wind_size);
        fogd_spectrum = conv(gaussFilter, enhanced_spec_peaks);
        fogd_spectrum = smooth(fogd_spectrum);

        % Apply smoothing several times x20
        % Eliminate Spurious VOPs
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

        % Normalize the FOGD Spectrum
        y = normalize(real(fogd_spectrum));
        y = y - min(y);
        y = y / max(y);

        % Create a time axis for FOGD Spectrum
        time_fogd = linspace(0, num_frames*frame_size*frame_over, length(fogd_spectrum));
        x = time_fogd;

        % Vowel Onset Points
        [pos_peaks, pos_peaks_idx] = findpeaks(y);

        % If only one vowel onset point exists
        if (length(pos_peaks_idx) >= 2)
            start_idx = floor(x(pos_peaks_idx(1))*fs);
            end_idx = floor(x(pos_peaks_idx(length(pos_peaks_idx)))*fs);

            l = 1;
            vowel_region= [];
            for k=start_idx:end_idx
                vowel_region(l) = data(k);
                l = l + 1;
            end

            audiowrite(new_file, vowel_region, fs);

        else
            audiowrite(new_file, data, fs);
        end
    end
end
