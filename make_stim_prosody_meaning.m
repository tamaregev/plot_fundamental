%% Definitions

dropbox_dir = '/Users/tamaregev/Dropbox/';
exp_dir = [dropbox_dir 'postdoc/Fedorenko/Prosody/Prosody-meaning'];
stim_folder = [exp_dir '/stimuli/'];
GH_folder = [exp_dir '/GitHub/prosody_meaning'];
cd(GH_folder)
CommonResources = [dropbox_dir 'MATLAB/lab/CommonResources'];
myFunctions = '/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions';
addpath(genpath(myFunctions))
addpath(CommonResources)
addpath(genpath(GH_folder));%from the github repository - includes STRAIGHT
speakers = {'WP','HK','HS'};
sentences = {[4 9 10 15 19 20],[1 2 3 5 6 7 8],[11 12 13 14 16 17 18]};%per speaker
%% Trim
ispeaker = 1;
speaker = speakers{ispeaker};
figure
for sent=sentences{ispeaker}
%for sent = 8
    raw_folder = [stim_folder speaker '/raw/Sent' num2str(sent) '/'];
    trim_folder = [stim_folder speaker '/trim/Sent' num2str(sent) '/'];
    
    if ~exist(trim_folder,'dir')
        mkdir(trim_folder)
    end
    
    Files = dir([raw_folder, '/*.m4a']);

    for ifile = 1:length(Files)
        filename = Files(ifile).name;
        [y,fs]=audioread([raw_folder filename]);
        if strcmp(speaker,'HK')
            y=mean(y,2);
        end
        soundsc(y,fs)
        %trim
        plot(y)
        [x,~]=ginput(2);
        x0=round(x(1));xend=round(x(2));
        fprintf(['x0=' num2str(x0) ' xend=' num2str(xend) '\n'])
        y = y(x0:xend);

        %add  50 ms fade in and out;
        win = 0.05;%50 ms
        nsamps = round(fs*win);
        y(1:nsamps) = y(1:nsamps).*linspace(0,1,nsamps)';
        y(end-nsamps+1:end) = y(end-nsamps+1:end).*linspace(1,0,nsamps)';
        
        plot(y)
        soundsc(y,fs)
        pause(3)
        audiowrite([trim_folder filename(1:end-4) '.wav'],y,fs)
    end
end
%% fix
ispeaker = 3;
speaker = speakers{ispeaker};
for sent=sentences{ispeaker}
    trim_folder = [stim_folder speaker '/trim/Sent' num2str(sent) '/'];
        
    Files = dir([trim_folder, '/*.m4a.wav']);

    for ifile = 1:length(Files)
        file = [trim_folder Files(ifile).name];
        [y,fs]=audioread([file]);
        audiowrite([file(1:end-8) '.wav'],y,fs)
        delete(file)
    end
end
%% Plot f0
params.th_f0score = 0.8;
params.th_df = 95;%Hz
params.conv = 5;

[h,f0_in,t0_in] = plotStimF0(y,fs,params, true);
hst=suptitle(strrep(filename,'_',' '));
set(hst,'fontsize',16)
%% Synth + manipulate
% https://docs.google.com/document/d/16xgOOVTrYqTQIVgNmAheyNmt1VvdPhlqmS4GwwnkkU4/edit
ispeaker = 1;
speaker = speakers{ispeaker};
for sent = sentences{ispeaker} 
    synth_folder = [stim_folder speaker '/synth/Sent' num2str(sent) '/'];
    if ~exist(synth_folder,'dir')
        mkdir(synth_folder)
    end   
end
sent = 10;
trim_folder = [stim_folder speaker '/trim/Sent' num2str(sent) '/'];
cd(trim_folder)
TandemSTRAIGHThandler
%I saved the synth versions with the prefix 'Syn' (this was automatically
%assigned by STRAIGHT) into the synth_folder manually from STRAIGHT GUI
%% Match length
speaker = speakers{1};
sent = 20;
synth_folder = [stim_folder speaker '/synth/Sent' num2str(sent) '/'];
category = 'N';version = 1;
prefix = 'Syn';
filename = [prefix num2str(sent) '_' category '_' num2str(version)];
[y,fs]=audioread([synth_folder filename '.wav']);
nSamps = length(y);

suffix = '_Q';
filename = [prefix num2str(sent) '_' category '_' num2str(version) suffix];
[y,fs]=audioread([synth_folder filename '.wav']);
y=y(1:nSamps);
audiowrite([synth_folder filename '.wav'],y,fs)
%% Morph
% https://memcauliffe.com/straight_workshop/morphing.html
speaker = 'WP';
sent = 20;
synth_folder = [stim_folder speaker '/synth/Sent' num2str(sent) '/'];
cd(synth_folder)
morph_folder = [stim_folder speaker '/morph/Sent' num2str(sent) '/'];
if ~exist(morph_folder,'dir')
    mkdir(morph_folder)
end
MorphingMenu
%% plot morph continuum
speaker = 'WP';
sent = 20;
pair1 = 'N';pair2 = 'Q';
name = [num2str(sent) '_' pair1 '_' pair2 ];
nmorphs = 4;
folder = [stim_folder speaker '/morph/Sent' num2str(sent) '/' name '/' name '_' num2str(nmorphs) ];

params.th_f0score = 0.8;
params.th_df = 95;%Hz
params.conv = 5;
plotf0morphs(params, folder, name);
legend({'morph #1','morph #2','morph #3','morph #4'},'location','se')
%% Test morphing continuum
% here I synthesized harmonic tones iwht a single f0, and morphed them
% with STRAIGHT to discover the spacing of the morph and the logarithmic
% spacing
f1=180;
f2=360;
dur = 2;ramp = [50 50];silence = [0 0];reduction = [0.5];fs=48000;
harmonics = [1:10];
phases = zeros(size(harmonics));
energy = reduction;

% I uploaded SynthHarmonicTone and SynthPureTone to the github in case
% anyone wants to play with them. STRAIGHT was able to detect well the
% pitch of the harmonic but not the pure tones.
[ data1, t1 ] = SynthHarmonicTone( f1, harmonics, phases, dur, ramp, silence, energy, fs );
[ data2, t2 ] = SynthHarmonicTone( f2, harmonics, phases, dur, ramp, silence, energy, fs );
audiowrite(['Harmonic' num2str(f1) '.wav'],data1,fs)
audiowrite(['Harmonic' num2str(f2) '.wav'],data2,fs)

%[ data1, t1 ] = SynthPureTone( f1, dur, ramp, silence, reduction, fs );
%[ data2, t2 ] = SynthPureTone( f2, dur, ramp, silence, reduction, fs );
%audiowrite(['sin' num2str(f1) '.wav'],data1,fs)
%audiowrite(['sin' num2str(f2) '.wav'],data2,fs)

plot(t1,data1)
hold on
plot(t2,data2)

MorphingMenu
% I morphed 100  
folder = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/morphSin/HarmonicContinuum5';
[h,f0_in,t0_in] = plotf0morphs(params, folder, 'Harmonic Continuum 5');
