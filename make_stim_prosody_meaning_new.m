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
sentences = 1:27;
%% Trim

h=figure;
set(h,'units','normalized','outerposition',[0 0 1 1])
%for sent = 8

    raw_folder = [stim_folder '/raw/'];
    trim_folder = [stim_folder '/trim/'];
    
    if ~exist(trim_folder,'dir')
        mkdir(trim_folder)
    end
    
    Files = dir([raw_folder, '/*.m4a']);

    for ifile = 1:length(Files)
        filename = Files(ifile).name;
        [y,fs]=audioread([raw_folder filename]);
        y=mean(y,2);
        soundsc(y,fs)
        clf
        plot(y)
        title(filename)

        [x,~]=ginput(2);
        x0=round(x(1));xend=round(x(2));
        fprintf(['x0=' num2str(x0) ' xend=' num2str(xend) '\n'])
        %trim
        y = y(x0:xend);

        %add  50 ms fade in and out;
        ramp_dur = 50;%ms
        y=hann_fade(y,ramp_dur,fs);
        
        plot(y)
        soundsc(y,fs)
        pause(3)
        audiowrite([trim_folder filename(1:end-4) '.wav'],y,fs)
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

for sent = sentences
    synth_folder = [stim_folder '/synth_raw/Sent' num2str(sent) '/'];
    if ~exist(synth_folder,'dir')
        mkdir(synth_folder)
    end   
end
trim_folder = [stim_folder '/trim/'];
cd(trim_folder)
TandemSTRAIGHThandler
%Saved the synth versions with the prefix 'Syn' (this was automatically
%assigned by STRAIGHT) into the synth_raw folder manually from STRAIGHT GUI
%% Match length
% and ramp up and down (fade in fade out) with a hanning window
%saved the synth folder before this section as: synth_raw
ramp_dur = 50; % ms

for sent = sentences
    synth_folder = [stim_folder '/synth_match/Sent' num2str(sent) '/'];
    if ~exist(synth_folder,'dir')
        mkdir(synth_folder)
    end   
end
for sent=sentences
    synth_in_folder = [stim_folder '/synth_raw/Sent' num2str(sent) '/'];
    %find length of shortest sound file:
    Files = dir([synth_in_folder, '/*.wav']);
    for i=1:length(Files)
        if i==1
            nsamps=nan(size(Files));
        end
        [y,fs]=audioread([synth_in_folder Files(i).name]);    
        nsamps(i) = length(y);
    end
    nsamp = min(nsamps);
    disp(sent)
    disp({Files(:).name})
    disp(nsamps')
    
    synth_out_folder = [stim_folder '/synth_match/Sent' num2str(sent) '/'];
    %load all other sound files, make same length and save to synth match folder:
    Files = dir([synth_folder, '/*.wav']);
    for i=1:length(Files)
        [y,fs]=audioread([synth_in_folder Files(i).name]);    
        y=y(1:nsamp);
        y=hann_fade(y,ramp_dur,fs);
        audiowrite([synth_out_folder Files(i).name],y,fs)
    end
    
    Files = dir([synth_out_folder, '/*.wav']);
    for i=1:length(Files)
        if i==1
            nsamps=nan(size(Files));
        end
        [y,fs]=audioread([synth_out_folder Files(i).name]);    
        nsamps(i) = length(y);
    end
    nsamp = min(nsamps);
    disp(sent)
    disp({Files(:).name})
    disp(nsamps')
    
    
end
%% make synth folder
for sent = sentences
    synth_folder = [stim_folder '/synth/Sent' num2str(sent) '/'];
    if ~exist(synth_folder,'dir')
        mkdir(synth_folder)
    end   
end
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
