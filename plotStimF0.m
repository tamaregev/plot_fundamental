function [h,f0_in,t0_in] = plotStimF0(y,fs,params, showProcess)
% plotStimF0 Calculates fundamental frequency using STRAIGHT, cleans it and plots signal+f0+spectrogram
% 
% Uses: STRAIGHT (Hideki Kawahara)
%       myspectrogram (by Kamil Wojcicki, Mathworks file exchange) - https://www.mathworks.com/matlabcentral/fileexchange/29596-speech-spectrogram
%
% Inputs:  y and fs are the outpts of [y,fs] = audioread('<soundFile>');
%
%          params is a structure
%          Example:
%               params.th_f0score = 0.75;%threshold f0 score for plotting f0 (score is 0 to 1)
%               params.th_df = 95;%Percentile maximal f0 jump for plotting f0
%               params.conv = 5;
%
%           showProcess is an optional logical varible. If showProcess = true a graph showing the f0 cleaning process will be displayed. Default: true. 


%
% June 2020, Tamar Regev, initial version
% June 15 2020, Tamar Regev
%                - Changed f0 cleaning process.
%                    input params subfield params.th_f0power replaced by params.th_f0score 
%                - Added plot to show process of f0 cleaning + optional
%                input variable showProcess
% 
%% Definitions
function_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/GitHub/plot_fundamental';
STRAIGHTFolder = [function_dir '/STRAIGHT/TandemSTRAIGHTmonolithicPackage010'];
addpath(STRAIGHTFolder)
addpath([function_dir '/myspectrogram'])

if ~exist('showProcess','var')
    showProcess=true;
end
    
%% Analysis
t = (0:length(y)-1)/fs;

r = exF0candidatesTSTRAIGHTGB(y,fs); % Extract F0 information
rc = autoF0Tracking(r,y); % Clean F0 trajectory by tracking
rc.vuv = refineVoicingDecision(y,rc);
q = aperiodicityRatioSigmoid(y,rc,1,2,0); % aperiodicity extraction
%f = exSpectrumTSTRAIGHTGB(y,fs,q); 

%f0power = q.f0candidatesPowerMap(1,:)';
f0score = q.f0CandidatesScoreMap(1,:)';
f0_in = rc.f0;
t0_in=rc.temporalPositions;

%select areas with high f0 power   
% th_f0power=max(q.f0candidatesPowerMap(1,:)')/params.th_f0power;
%th_f0power=prctile(q.f0candidatesPowerMap(1,:),45);

%select areas with high f0 score   
th_f0score =  params.th_f0score;%typically 0.75


%select areas with big jumps:
th_df=prctile(diff(f0_in),params.th_df);
df0=nan(size(f0_in));
df0(1:end-1)=abs(diff(f0_in));

if showProcess
    figure;
    subplot(2,1,1)
    plot(rc.temporalPositions,f0score,'.');grid on;hold on
    plot([rc.temporalPositions(1),rc.temporalPositions(end)],[th_f0score,th_f0score])
    xlim([0 t(end)])
    title('F0 score + threshold')
    set(gca,'fontsize',14)
    
    subplot(2,1,2)
    plot(t0_in,f0_in,'.','markersize',14);grid on;hold on
    ylabel('f0 (Hz)');
    ylim([50 350])
    xlim([0 t(end)])
    plot(t0_in(df0>th_df),f0_in(df0>th_df),'.r','markersize',18)
    plot(t0_in(f0score<th_f0score),f0_in(f0score<th_f0score),'.g')
    title('raw + cleaned F0')
    legend('clean F0','rejected due to large df','rejected due to low score','fontsize',12,'location','sw')
    set(gca,'fontsize',14)
    
end

%smooth:
convolution = ones(1,params.conv)/params.conv;
f0_in=conv(f0_in,convolution,'same');

%clean f0:
f0_in(df0>th_df)=nan;
t0_in(df0>th_df)=nan;
f0_in(f0score<th_f0score)=nan;
t0_in(f0score<th_f0score)=nan;
    

%% plot
    h=figure;
    
    %%% plot signal
    subplot(3,1,1); plot(t,y);xlim([0 t(end)]); grid on
    set(gca,'fontsize',14)
        
    %%% plot fundamental 
    subplot(3,1,2)
    %plot(t0_in,f0_in,'.');grid on
    semilogy(t0_in,f0_in,'.');grid on
    ylabel('f0 (Hz)');
    ylim([100 350])
    xlim([0 t(end)])
    set(gca,'fontsize',14)
    
    %%% plot spectrogram
    subplot(3,1,3)
    [hs, T, F] = myspectrogram(y, fs, [22 1], @hamming, 2048, [-59 -1], false, 'default', false, 'per'); % or be quite specific about what you want
    ylim([0 2000])
    xlim([0 t(end)])
    colormap('parula')
    xlabel('Time (s)')
    set(gca,'fontsize',14)

end

