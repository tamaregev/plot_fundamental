%plot_results_acoustic_analysis_morphs
function plot_results_acoustic_analysis_morphs(morphFolder,delta_pitch_tag)
    function_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/GitHub/prosody_meaning';
    addpath(genpath(function_dir))
    T = readtable([morphFolder filesep 'analysis_dp_' delta_pitch_tag '.csv']);
    [ss,I] = sort(T.sent);
    Ts = T(I,:);
    %% Setting variables - critical sentences
    critSent = 1:18;
    NcritSent = numel(critSent);
    measures = {'maxpitch','sumpitch','maxpitch_O2','sumpitch_O2','maxloud','sumloud','maxloud_O2','sumloud_O2'};
    conditions = {'No','F','Q'};
    stepmorphs = 2:4;
    Y=nan(numel(stepmorphs),numel(conditions),numel(measures),NcritSent);
    for meas = 1:numel(measures)
        for con = 1:length(conditions)
            for morph = 1:length(stepmorphs)
                Y(morph,con,meas,:) = Ts.(measures{meas})(strcmp(Ts.condition,conditions{con}) & ismember(Ts.stepmorph,stepmorphs(morph)) & ismember(Ts.sent,critSent));
            end
        end
    end
    Y_mean = nanmean(Y,4);
    Y_stde = nanstd(Y,0,4)/sqrt(NcritSent);

    %%
    addpath('/Users/tamaregev/Dropbox/MATLAB/lab/CommonResources');

    ERPfigure
    switch delta_pitch_tag
        case 'difference'
            yUnits = {'Hz','Hz','Hz/s','Hz/s','dB','dB','dB/s','dB/s'};
        case 'fraction'
            yUnits = {'% (dHz/Hz)','% (dHz/Hz)','Hz/s','Hz/s','dB','dB','dB/s','dB/s'};
    end
    for meas = 1:length(measures)
        subplot(2,4,meas)
        h{meas} = bar(Y_mean(:,:,meas)');
        title(strrep(measures{meas},'_',' '),'fontsize',18)
        set(gca, 'XTickLabel', conditions);
        ylabel(yUnits{meas})
        for ii=1:numel(h{meas})%goes by color = morph
            xlocs = h{meas}(ii).XEndPoints;
            RGB = h{meas}(ii).FaceColor;
            for iii=1:numel(xlocs)%cond
                hold on
                scatter(repmat(xlocs(iii),[NcritSent,1]),squeeze(Y(ii,iii,meas,:)),30,RGB)
                ts = textscatter(repmat(xlocs(iii),[NcritSent,1]),squeeze(Y(ii,iii,meas,:)),string(critSent),'TextDensityPercentage',100);
            end
        end
        set(gca,'fontsize',18)
    end
    legend({'morph 1-2','morph 1-3','morph 1-4'})
    suptitle({'Acoustic differences',morphFolder})