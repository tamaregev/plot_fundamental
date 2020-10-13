%plot_results_acoustic_analysis_morphs
morphFolder = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/stimuli/morph';
T=readtable([morphFolder filesep 'analysis.csv']);

%% Setting variables - critical sentences
critSent = 1:18;
NcritSent = numel(critSent);
measures = {'maxpitch','sumpitch','maxloud','sumloud'};
conditions = {'No','F','Q'};
stepmorphs = 2:4;
Y=nan(numel(stepmorphs),numel(conditions),numel(measures),NcritSent);
for meas = 1:numel(measures)
    for con = 1:length(conditions)
        for morph = 1:length(stepmorphs)
            Y(morph,con,meas,:) = T.(measures{meas})(strcmp(T.condition,conditions{con}) & ismember(T.stepmorph,stepmorphs(morph)) & ismember(T.sent,critSent));
        end
    end
end
Y_mean = nanmean(Y,4);
Y_stde = nanstd(Y,0,4)/sqrt(NcritSent);

%%
addpath('/Users/tamaregev/Dropbox/MATLAB/lab/CommonResources');

ERPfigure
yUnits = {'Hz','Hz','dB','dB'};
for meas = 1:4
    subplot(2,2,meas)
    h{meas} = bar(Y_mean(:,:,meas)');
    title(measures{meas},'fontsize',18)
    set(gca, 'XTickLabel', conditions);
    ylabel(yUnits{meas})
    for ii=1:numel(h{meas})%goes by color = morph
        xlocs = h{meas}(ii).XEndPoints;
        RGB = h{meas}(ii).FaceColor;
        for iii=1:numel(xlocs)%cond
            hold on
         %   scatter(repmat(xlocs(iii),[NcritSent,1]),squeeze(Y(ii,iii,meas,:)),30,RGB)
            ts = textscatter(repmat(xlocs(iii),[NcritSent,1]),squeeze(Y(ii,iii,meas,:)),string(critSent),'TextDensityPercentage',100);
        end
    end
    set(gca,'fontsize',14)
end
legend({'morph 1-2','morph 1-3','morph 1-4'})
