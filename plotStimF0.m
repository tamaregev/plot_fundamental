function [h,f0_in,t0_in] = plotStimF0(y,fs,params)
% plotStimF0 Calculates fundamental frequency using STRAIGHT, cleans it and plots signal+f0+spectrogram
% 
% Uses: STRAIGHT (Hideki Kawahara)
%       myspectrogram (by Kamil Wojcicki, Mathworks file exchange) - https://www.mathworks.com/matlabcentral/fileexchange/29596-speech-spectrogram
%
% Inputs: params is a structure
%          Example:
%               params.th_f0power = 0.017;threshold f0 power for plotting f0
%               params.th_df = 10;%Hz maximal f0 jump for plotting f0
%               params.conv = 5;

%
% Tamar Regev, June 2020
%% Definitions
dropbox_dir = '/Users/tamaregev/Dropbox/';
STRAIGHTFolder = [dropbox_dir 'MATLAB/lab/STRAIGHT/STRAIGHT/2015/TandemSTRAIGHTmonolithicPackage010'];
addpath(STRAIGHTFolder)
addpath([dropbox_dir filesep 'MATLAB/lab/myFunctions/downloaded/myspectrogram'])

%% Analysis
t = (0:length(y)-1)/fs;

r = exF0candidatesTSTRAIGHTGB(y,fs); % Extract F0 information
rc = autoF0Tracking(r,y); % Clean F0 trajectory by tracking
rc.vuv = refineVoicingDecision(y,rc);
q = aperiodicityRatioSigmoid(y,rc,1,2,0); % aperiodicity extraction
f = exSpectrumTSTRAIGHTGB(y,fs,q); 

    f0power = q.f0candidatesPowerMap(1,:)';
    f0_in = rc.f0;
    t0_in=rc.temporalPositions;

    %select areas with high f0 power
    th_f0power=params.th_f0power;
    f0_in(f0power<th_f0power)=nan;
    t0_in(f0power<th_f0power)=nan;
    
    %clean f0:
    %exclude big jumps:
    th_df = params.th_df;%Hz
    df0=nan(size(f0_in));
    df0(1:end-1)=abs(diff(f0_in));
    f0_in(df0>th_df)=nan;
    t0_in(df0>th_df)=nan;
    %exclude isolated points:
%     
    %smooth:
    convolution = ones(1,params.conv)/params.conv;
   f0_in=conv(f0_in,convolution,'same');
   

%% plot
h=figure;
     subplot(3,1,1); plot(t,y);xlim([0 t(end)]); grid on
     set(gca,'fontsize',14)
%      subplot(7,1,2); plot(t0,f0); grid on
%         ylabel('Pitch (Hz)')
%         ylim([50 250])
%         xlim([0 t(end)])
        
        subplot(3,1,2)
        %%% plot fundamental 
%         plot(rc.temporalPositions,rc.f0);grid on
%         set(gca,'fontsize',14);
%         xlabel('Time (s)')
%         ylabel('fundamental frequency (Hz)');
%         ylim([50 250])
%         xlim([0 t(end)])

        plot(t0_in,f0_in,'linewidth',2);grid on
        
        ylabel('f0 (Hz)');
        ylim([50 250])
        xlim([0 t(end)])
        set(gca,'fontsize',14)
%        subplot(7,1,4);hold off
%         %figure
%         plot(rc.temporalPositions,q.f0candidatesPowerMap(1,:)') 
%         xlim([0 t(end)])
%         grid on
%         
%         subplot(7,1,5);hold off
%         plot(rc.temporalPositions,q.periodicityLevel) 
%         xlim([0 t(end)])
%         grid on

        subplot(3,1,3)
        [hs, T, F] = myspectrogram(y, fs, [22 1], @hamming, 2048, [-59 -1], false, 'default', false, 'per'); % or be quite specific about what you want
        ylim([0 2000])
        xlim([0 t(end)])
        colormap('parula')
        xlabel('Time (s)')
        set(gca,'fontsize',14)
%         subplot(7,1,7)
%         ss=sum(hs.CData);
%         plot(T,ss)
%         xlim([T(1) T(end)])
%         grid on


end

