%% Definitions

dropbox_dir = '/Users/tamaregev/Dropbox/';
stim_dir = [dropbox_dir 'postdoc/Fedorenko/Prosody/Prosody-meaning/stimuli/'];
CommonResources = [dropbox_dir 'MATLAB/lab/CommonResources'];
myFunctions = '/Users/tamaregev/Dropbox/MATLAB/lab/myFunctions';
addpath(genpath(myFunctions))
addpath(CommonResources)
addpath(genpath('GitHub/prosody_meaning'));%from the github repository - includes STRAIGHT
speakers = {'WP','HK'};

%% Trim
speaker = speakers{1};
sent = 20;
raw_folder = [stim_dir speaker '/raw/Sent' num2str(sent) '/'];
trim_folder = [stim_dir speaker '/trim/Sent' num2str(sent) '/'];
if ~exist(trim_folder,'dir')
    mkdir(trim_folder)
end

category = 'F2';version = 2;
filename = [num2str(sent) '_' category '_' num2str(version)];
[y,fs]=audioread([raw_folder filename '.m4a']);
soundsc(y,fs)

%trim
%plot(y)
%x0=63420; xend=192800;4 N 1
%x0=45320;xend=156500;
x0=48700;xend=147200;

y = y(x0:xend);

%add  50 ms fade in and out;
win = 0.05;%50 ms
nsamps = round(fs*win);
y(1:nsamps) = y(1:nsamps).*linspace(0,1,nsamps)';
y(end-nsamps+1:end) = y(end-nsamps+1:end).*linspace(1,0,nsamps)';

plot(y)
soundsc(y,fs)

audiowrite([trim_folder filename '.wav'],y,fs)
%% Analyze f0
params.th_f0score = 0.8;
params.th_df = 95;%Hz
params.conv = 5;

[h,f0_in,t0_in] = plotStimF0(y,fs,params, true);
hst=suptitle(strrep(filename,'_',' '));
set(hst,'fontsize',16)
%% Synth + manipulate
% https://docs.google.com/document/d/16xgOOVTrYqTQIVgNmAheyNmt1VvdPhlqmS4GwwnkkU4/edit
TandemSTRAIGHThandler
%I saved the synth versions with the prefix 'Syn' (this was automatically
%assigned by STRAIGHT)
%% Match length
speaker = speakers{1};
sent = 20;
synth_folder = [stim_dir speaker '/synth/Sent' num2str(sent) '/'];
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
MorphingMenu
%% plot morph continuum

folder = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/morphSin/HarmonicContinuum5';
[h,f0_in,t0_in] = plotf0morphs(params, folder, 'Harmonic Continuum 5');

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
