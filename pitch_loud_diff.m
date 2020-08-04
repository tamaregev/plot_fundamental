function [delta_pitch, delta_loud] = pitch_loud_diff(stimulus_1,stimulus_2,params)
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
% results_f0 = [];
% results_db = [0 0];

for ii = 1:2
[y,fs] = audioread(stimuli(ii));
% t = (0:length(y)-1)/fs; % for graphing
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
        results_f0 = [results_f0,f0_in];
    end
%%
% wav_Pa = y * params.my_cal_factor;
% smooth_sec = 0.125;  %"FAST" SPL is 1/8th of second.  "SLOW" is 1 second;
% smooth_Hz = 1/smooth_sec;
% 
% [b,a]=butter(1,smooth_Hz/(fs/2),'low');  %design a Low-pass filter
% wav_env_Pa = sqrt(filter(b,a,wav_Pa.^2));  %rectify, by squaring, and low-pass filter
% 
% %compute SPL
% Pa_ref = 20e-6;  %reference pressure for SPL in Air
% SPL_dB = 10.0*log10( (wav_env_Pa ./ Pa_ref).^2 ); % 10*log10 because signal is squared
%     if ii == 1
%         results_db(ii) = SPL_dB;
%         results_db(ii+1) = t0_in;
%     else    
%         results_db(ii+1) = SPL_dB;
%         results_db(ii+2) = t0_in;
%     end  
end
delta_pitch = results_f0(:,1)./results_f0(:,2);
delta_loud = 1; %placeholder
% delta_loud = results_db(:,1) / results_db(:,2);
end