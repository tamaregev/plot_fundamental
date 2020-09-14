function [h,f0_in,t0_in] = plotf0morphs(params, folder, tag)
% plotf0morphs calculates fundamental frequency using STRAIGHT, cleans it and
% plots f0 of all wav files of continuum
% 
% Uses: STRAIGHT (Hideki Kawahara)
%       
% Inputs:  params is a structure
%          Example:
%               params.th_f0score = 0.75;threshold f0 score for plotting f0 (score is 0 to 1)
%               params.th_df = 95;%Percentile maximal f0 jump for plotting f0
%               params.conv = 5;
%         
%          folder is a string, the full path to the folder
%          containing stimuli.
%

%
% Based on plotStimF0 (Tamar Regev)
% June 22, 2020 Julie Meng
% July 8, 2020 Tamar Regev: changed inputs such that folder is the full
% folder and tag is the title for the figure.
% Aug 11 2020, Tamar Regev
%               - Fixed bug added abs in th_df=prctile(abs(diff(f0_in)),params.th_df);

%% Definitions
function_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/GitHub/prosody_meaning';
addpath(genpath(function_dir))

%% Analysis
Files = dir([folder, '/*.wav']);

FileNames = string.empty;
for k=1:length(Files)
    %if strfind(Files(k).name, 'wav')
   FileNames = [FileNames Files(k).name];
end

for k=1:length(FileNames)
    disp(FileNames(k))
end

for index = 1:length(FileNames)
    [y,fs] = audioread([folder filesep FileNames{index}]); %read audio file
    
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
    th_f0score = params.th_f0score;%typically 0.75

    %select areas with big jumps:
    th_df=prctile(abs(diff(f0_in)),params.th_df);
    df0=nan(size(f0_in));
    df0(1:end-1)=abs(diff(f0_in));

    %smooth:
    convolution = ones(1,params.conv)/params.conv;
    f0_in=conv(f0_in,convolution,'same');

    %clean f0:
    f0_in(df0>th_df)=nan;
    t0_in(df0>th_df)=nan;
    f0_in(f0score<th_f0score)=nan;
    t0_in(f0score<th_f0score)=nan;
    
    if index == 1
        h=figure;
        
        %%% plot fundamental 
        semilogy(t0_in,f0_in,'.');grid on

        hold on
     
    elseif index == length(FileNames)
        semilogy(t0_in,f0_in,'.')
        ylabel('f0 (Hz)');
        xlabel('Time (s)')
        ylim([100 450])
        xlim([0 t(end)])
        title(strrep(tag, "_", " "))
        set(gca,'fontsize',14)
        hold off
        
    else
        plot(t0_in,f0_in,'.')
    end
        
    
end

end