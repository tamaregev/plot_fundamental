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
    Files = dir([synth_in_folder, '/*.wav']);
    for i=1:length(Files)
        [y,fs]=audioread([synth_in_folder Files(i).name]);    
        y=y(1:nsamp);
        y=hann_fade(y,ramp_dur,fs);
        audiowrite([synth_out_folder Files(i).name],y,fs)
    end
    
end
%check and match that stimuli in synth match folder are the same length
for sent=sentences
    synth_out_folder = [stim_folder '/synth_match/Sent' num2str(sent) '/'];

    Files = dir([synth_out_folder, '/*.wav']);
    for i=1:length(Files)
        if i==1
            nsamps=nan(size(Files));
        end
        [y,fs]=audioread([synth_out_folder Files(i).name]);    
        nsamps(i) = length(y);
    end
    nsamp = min(nsamps);
    disp([num2str(sent) ' before'])
    disp({Files(:).name})
    disp(nsamps')
    for i=1:length(Files)
        [y,fs]=audioread([synth_out_folder Files(i).name]);    
        y=y(1:nsamp);        
        nsamps(i) = length(y);
        y=hann_fade(y,ramp_dur,fs);
        audiowrite([synth_out_folder Files(i).name],y,fs)
    end
    disp([num2str(sent) ' after'])
    disp({Files(:).name})
    disp(nsamps')
end

%check that stimuli in synth folder are the same length
for sent=sentences
    synth_folder = [stim_folder '/synth/Sent' num2str(sent) '/'];
    %find length of shortest sound file:
    Files = dir([synth_folder, '/*.wav']);
    for i=1:length(Files)
        if i==1
            nsamps=nan(size(Files));
        end
        [y,fs]=audioread([synth_folder Files(i).name]);    
        nsamps(i) = length(y);
    end
    nsamp = min(nsamps);
    disp(sent)
    disp({Files(:).name})
    disp(nsamps')
    
    for i=1:length(Files)
        [y,fs]=audioread([synth_folder Files(i).name]);    
        y=y(1:nsamp);        
        nsamps(i) = length(y);
        y=hann_fade(y,ramp_dur,fs);
        audiowrite([synth_folder Files(i).name],y,fs)
    end
    disp([num2str(sent) ' after'])
    disp({Files(:).name})
    disp(nsamps')
end
%% make synth folder
for sent = sentences
    synth_folder = [stim_folder '/synth/Sent' num2str(sent) '/'];
    if ~exist(  synth_folder,'dir')
        mkdir(synth_folder)
    end   
end
%% Morph
% https://memcauliffe.com/straight_workshop/morphing.html
sent = 2;
synth_folder = [stim_folder '/synth/Sent' num2str(sent) '/'];
cd(synth_folder)
morph_folder = [stim_folder '/morph/Sent' num2str(sent) '/'];
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
%% Make masks - avg
mask_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/stimuli/mask_avg';
nMasks = 10;%how many different masks to make
nStim = 50;%per mask
[mask_all] = makeMasks_avg(mask_dir,nMasks,nStim);

%% Make masks - pink
mask_dir = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/stimuli/mask_pink';
pinkNoise = dsp.ColoredNoise(1,44.1e3,1);
rng default;
x = pinkNoise();
soundsc(x,fs)
audiowrite([mask_dir filesep 'pink_noise.wav'],x,fs)

%% Fix new synth (ballance acoustics)
sent=18;
cond = 'Q';
ramp_dur = 50; % ms
suffix = '_balp';

filename_new = ['Syn' num2str(sent) '_' cond suffix '.wav'];
filename_old = ['Syn' num2str(sent) '_' cond '.wav'];
synth_raw_folder = [stim_folder '/synth_raw/Sent' num2str(sent) '/'];
synth_match_folder = [stim_folder '/synth_match/Sent' num2str(sent) '/'];
synth_folder = [stim_folder '/synth/Sent' num2str(sent) '/'];

cd(synth_raw_folder)

%copy to synth_match
copyfile([synth_raw_folder filename_new],[synth_match_folder filename_new]);
cd(synth_match_folder)

%check and match that stimuli in synth match folder are the same length
[y,fs]=audioread([synth_match_folder filename_new]);    
nsamps = length(y);
disp(['new file length = ' num2str(nsamps) ' samps'])

[y,fs]=audioread([synth_match_folder filename_old]);    
nsamps = length(y);
disp(['old file length = ' num2str(nsamps) ' samps'])

%calc length of all files in matched folder
Files = dir([synth_match_folder, '/*.wav']);
for i=1:length(Files)
    if i==1
        nsamps=nan(size(Files));
    end
    [y,fs]=audioread([synth_match_folder Files(i).name]);    
    nsamps(i) = length(y);
end
nsamp = min(nsamps);
disp({Files(:).name})
disp(nsamps')
    
%match their length , applying a hann window
for i=1:length(Files)
    [y,fs]=audioread([synth_match_folder Files(i).name]);    
    y=y(1:nsamp);        
    nsamps(i) = length(y);
    y=hann_fade(y,ramp_dur,fs);
    audiowrite([synth_match_folder Files(i).name],y,fs)
end

%verify that all have the same length
disp({Files(:).name})
disp(nsamps')

%copy into synth folder
cd(synth_folder)
copyfile([synth_match_folder filename_new],[synth_folder filename_new]);

%calc length of all files in synth folder
Files = dir([synth_folder, '/*.wav']);
for i=1:length(Files)
    if i==1
        nsamps=nan(size(Files));
    end
    [y,fs]=audioread([synth_folder Files(i).name]);    
    nsamps(i) = length(y);
end
disp({Files(:).name})
disp(nsamps')

%match their length , applying a hann window
for i=1:length(Files)
    [y,fs]=audioread([synth_folder Files(i).name]);    
    y=y(1:nsamp);        
    nsamps(i) = length(y);
    y=hann_fade(y,ramp_dur,fs);
    audiowrite([synth_folder Files(i).name],y,fs)
end

%verify that all have the same length
disp({Files(:).name})
disp(nsamps')

    %% calculate acoustics of new synth
filename_neutral = ['Syn' num2str(sent) '_N.wav'];

[y0,fs]=audioread([synth_folder filename_neutral]);    
[y1,fs]=audioread([synth_folder filename_old]);
[y2,fs]=audioread([synth_folder filename_new]);


params.th_f0score = 0.75; %threshold f0 power for plotting f0
params.th_df = 95; %Hz maximal f0 jump for plotting f0
params.conv = 5;
params.my_cal_factor = 1.0;  %the value for your system to convert the WAV into Pascals

%compare old to N
[delta_pitch, delta_loud, p1, p2, t_p, l1, l2, t_l] = pitch_loud_diff(y0,y1,fs,params,true);
suptitle('old')
dp1=diff(p1);dp2=diff(p2);
dl1=diff(l1);dl2=diff(l2);
                                            
sumpitch(1,1) = nansum(p2-p1);
maxpitch(1,1) = nanmax(p2-p1);
sumpitch_O2(1,1) = nansum(abs(dp2)-abs(dp1)); 
maxpitch_O2(1,1) = nanmax(abs(dp2)-abs(dp1));

sumloud(1,1) = nansum(delta_loud);
maxloud(1,1) = nanmax(delta_loud);
sumloud_O2(1,1) = nansum(abs(dl2)-abs(dl1)); 
maxloud_O2(1,1) = nanmax(abs(dl2)-abs(dl1));

%compare new to N
[delta_pitch, delta_loud, p1, p2, t_p, l1, l2, t_l] = pitch_loud_diff(y0,y2,fs,params,true);
suptitle('new')

dp1=diff(p1);dp2=diff(p2);
dl1=diff(l1);dl2=diff(l2);
                                            
sumpitch(2,1) = nansum(p2-p1);
maxpitch(2,1) = nanmax(p2-p1);
sumpitch_O2(2,1) = nansum(abs(dp2)-abs(dp1)); 
maxpitch_O2(2,1) = nanmax(abs(dp2)-abs(dp1));

sumloud(2,1) = nansum(delta_loud);
maxloud(2,1) = nanmax(delta_loud);
sumloud_O2(2,1) = nansum(abs(dl2)-abs(dl1)); 
maxloud_O2(2,1) = nanmax(abs(dl2)-abs(dl1));

condition = {[cond '_old'],[cond '_new']}';
stepmorph = [4,4]';
sentence = [sent sent]';

Tb = table(sentence,condition,stepmorph,sumpitch,maxpitch,sumpitch_O2,maxpitch_O2,sumloud,maxloud,sumloud_O2,maxloud_O2)


morphFolder = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/stimuli/morph';
T = readtable([morphFolder filesep 'analysis.csv']);
T(T.sent==sent & ismember(T.condition,'Q'),:)
%% test differences statistically
morphFolder = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/stimuli/morph_balanced';
T = readtable([morphFolder filesep 'analysis_dp_difference_.csv']);
maxpitch_Q=T(ismember(T.sent,[1:18]) & ismember(T.condition,'Q') & T.stepmorph==4,:).maxpitch;
maxpitch_No=T(ismember(T.sent,[1:18]) & ismember(T.condition,'No') & T.stepmorph==4,:).maxpitch;
maxloud_No=T(ismember(T.sent,[1:18]) & ismember(T.condition,'No') & T.stepmorph==4,:).maxloud;
maxloud_Q=T(ismember(T.sent,[1:18]) & ismember(T.condition,'Q') & T.stepmorph==4,:).maxloud;
sent_numbers = T(ismember(T.sent,[1:18]) & ismember(T.condition,'Q') & T.stepmorph==4,:).sent;

[H,P,CI,STATS] = ttest(maxpitch_Q,maxpitch_No)
[H,P,CI,STATS] = ttest(maxloud_Q,maxloud_No)

%% create a morph4to3 dir
source = [stim_folder '/morph_balanced'];
destination = [stim_folder '/morph_balanced_morph4to3'];

copyfile(source, destination)

D=dir(destination);
D = D(~ismember({D.name}, {'.', '..','.DS_Store','figures'})); % dir returns '.' and '..' (usually in first slot)
D = D([D.isdir]); %only look at subdirectories

conditions= {'Q','No','F'};

for ifolder = 1:length(D)
    for icond = 1:length(conditions)
        subfolder = [D(ifolder).name(5:end) '_' conditions{icond}];
        morph3 = [destination filesep D(ifolder).name filesep subfolder filesep subfolder '_003.wav']; 
        morph4 = [destination filesep D(ifolder).name filesep subfolder filesep subfolder '_004.wav'];        
        copyfile(morph3, morph4)
    end
end
%% listen reversed
soundsc(y0,fs)
soundsc(y1,fs)
soundsc(y2,fs)

soundsc(y0(end:-1:1),fs)
soundsc(y1(end:-1:1),fs)
soundsc(y2(end:-1:1),fs)
