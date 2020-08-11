function [delta_pitch, delta_loud, t, p1, p2, l1, l2] = pitch_loud_diff(stimulus_1,stimulus_2,params)
%               params.th_f0score = 0.75; %threshold f0 power for plotting f0
%               params.th_df = 95; %Hz maximal f0 jump for plotting f0
%               params.conv = 5;
%               params.my_cal_factor = 1.0;  %the value for your system to convert the WAV into Pascals

%%
function_dir = 'C:\Users\user\Downloads\prosody_meaning-master\prosody_meaning-master\';
STRAIGHTFolder = [function_dir '\STRAIGHT\TandemSTRAIGHTmonolithicPackage010'];
addpath(STRAIGHTFolder)
addpath([function_dir '\myspectrogram'])

%%
stim1 = string(stimulus_1);
stim2 = string(stimulus_2);
stimuli = [stim1, stim2];

for ii = 1:2
[y,fs] = audioread(stimuli(ii));
th_f0score =  params.th_f0score; %typically 0.75

r = exF0candidatesTSTRAIGHTGB(y,fs); % Extract F0 information
rc = autoF0Tracking(r,y); % Clean F0 trajectory by tracking
rc.vuv = refineVoicingDecision(y,rc);
q = aperiodicityRatioSigmoid(y,rc,1,2,0); % aperiodicity extraction

    f0score = q.f0CandidatesScoreMap(1,:)';
    f0_in = rc.f0;
    t0_in = rc.temporalPositions;

    th_df = prctile(diff(f0_in),params.th_df);
    df0 = nan(size(f0_in));
    df0(1:end-1) = abs(diff(f0_in));
    
    %smooth:
    convolution = ones(1,params.conv)/params.conv;
    f0_in = conv(f0_in,convolution,'same');
    
    %clean f0:
    f0_in(df0>th_df) = nan;
    t0_in(df0>th_df) = nan;
    f0_in(f0score<th_f0score) = nan;
    t0_in(f0score<th_f0score) = nan;
    if ii == 1
        results_f0 = [f0_in];
    else
        [dim1,~] = size(f0_in);
        resize = dim1 - 1;
        f0_in = f0_in(1:resize,:);
        results_f0 = [results_f0,f0_in];
    end
%%
wav_Pa = y * params.my_cal_factor;
smooth_sec = 0.125;  %"FAST" SPL is 1/8th of second.  "SLOW" is 1 second;
smooth_Hz = 1/smooth_sec;

[b,a]=butter(1,smooth_Hz/(fs/2),'low');  %design a Low-pass filter
wav_env_Pa = sqrt(filter(b,a,wav_Pa.^2));  %rectify, by squaring, and low-pass filter

%compute SPL
Pa_ref = 20e-6;  %reference pressure for SPL in Air
SPL_dB = 10.0*log10( (wav_env_Pa ./ Pa_ref).^2 ); % 10*log10 because signal is squared
    if ii == 1
        results_db = [SPL_dB];
        [init_dim,~] = size(SPL_dB);
    else
        [dim1,~] = size(SPL_dB);
        if init_dim > dim1
            baseSPL_dB = nan(init_dim,1);
            baseSPL_dB(1:dim1,1) = SPL_dB;
            results_db = [results_db,baseSPL_dB];
        elseif init_dim < dim1
            baseSPL_dB = nan(dim1,1);
            baseSPL_dB(1:init_dim,1) = results_db;
            results_db = [baseSPL_dB,SPL_dB];
        end
    end  
end
t = t0_in;
p1 = results_f0(:,1);
p2 = results_f0(:,2);

l1 = results_db(:,1);
l2 = results_db(:,2);

delta_pitch = results_f0(:,2)./results_f0(:,1);
delta_loud = results_db(:,2) - results_db(:,1);
end