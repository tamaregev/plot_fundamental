function [delta_pitch, delta_loud, p1, p2, t_p, l1, l2, t_l] = pitch_loud_diff(y1,y2,fs,params,plotFlag)

% INPUTS:
%         y1, y2 are column sound signals, must be same length
%         fs is the sample rate
%         params is a structure, example:
%               params.th_f0score = 0.75; %threshold f0 power for plotting f0
%               params.th_df = 95; %Hz maximal f0 jump for plotting f0
%               params.conv = 5;
%               params.my_cal_factor = 1.0;  %the value for your system to convert the WAV into Pascals
%         plotFlag is either true or false determines if function makes
%         some plots
%%
if ~exist('plotFlag','var')
    plotFlag = true;
end


% function_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/GitHub/prosody_meaning/';
% STRAIGHTFolder = [function_dir '/STRAIGHT/TandemSTRAIGHTmonolithicPackage010'];
% addpath(STRAIGHTFolder)
% addpath([function_dir '/myspectrogram'])


%%
stimuli = [y1, y2];

for ii = 1:2
    y = stimuli(:,ii);

    %% F0 trajectory
    th_f0score =  params.th_f0score; %typically 0.75

    r = exF0candidatesTSTRAIGHTGB(y,fs); % Extract F0 information
    rc = autoF0Tracking(r,y); % Clean F0 trajectory by tracking
    rc.vuv = refineVoicingDecision(y,rc);
    q = aperiodicityRatioSigmoid(y,rc,1,2,0); % aperiodicity extraction

    f0score = q.f0CandidatesScoreMap(1,:)';
    f0_in = rc.f0;
    t0_in = rc.temporalPositions;

    th_df = prctile(abs(diff(f0_in)),params.th_df);
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
        p1 = f0_in;
    elseif ii == 2
        p2 = f0_in;
    end
%% loudness trajectory
    wav_Pa = y * params.my_cal_factor;
    smooth_sec = 0.125;  %"FAST" SPL is 1/8th of second.  "SLOW" is 1 second;
    smooth_Hz = 1/smooth_sec;

    [b,a]=butter(1,smooth_Hz/(fs/2),'low');  %design a Low-pass filter
    wav_env_Pa = sqrt(filter(b,a,wav_Pa.^2));  %rectify, by squaring, and low-pass filter

    %compute SPL
    Pa_ref = 20e-6;  %reference pressure for SPL in Air
    SPL_dB = 10.0*log10( (wav_env_Pa ./ Pa_ref).^2 ); % 10*log10 because signal is squared

   
    if ii == 1
        l1 = SPL_dB;
    elseif ii == 2
        l2 = SPL_dB;
    end  
end

t_l = 0:1/fs:(length(l1)-1)/fs;t_l=t_l';
t_p = t0_in;
t_p_all = rc.temporalPositions;

%clean l1, l2 and t_l according to p1, p2 and t_p
nan_times = t_p_all(isnan(t_p));
tol = 0.01;
t_l(ismembertol(t_l,nan_times,tol)) = nan;
l1(isnan(t_l)) = nan;
l2(isnan(t_l)) = nan;

if plotFlag % plotting loudness trajectories
        figure;
        subplot(3,1,1);
        t_sec = ([1:size(y1)]-1)/fs;
        plot(t_sec,y1);
        xlabel('Time (sec)');
        ylabel('Pressure (Pa)');
        xlim([0 3])
        
        subplot(3,1,2)
        plot(t_sec,l1);
        xlabel('Time (sec)');
        ylabel('SPL (dB)');
        yl=ylim;ylim(yl(2)+[-80 0]);
        hold all
        plot(t_sec,l2);
        xlim([0 3])
        
        subplot(3,1,3)
        plot(t_p,p1);
        hold all
        plot(t_p,p2);
        xlim([0 3])
        
end

delta_pitch = p2 - p1;
delta_loud = l2 - l1;
end
