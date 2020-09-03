%% current folder should be main folder with subdirectories of synth files
params.th_f0score = 0.75; %threshold f0 power for plotting f0
params.th_df = 95; %Hz maximal f0 jump for plotting f0
params.conv = 5;
params.my_cal_factor = 1.0;  %the value for your system to convert the WAV into Pascals

D = dir; % A is a struct ... first elements are '.' and '..' used for navigation.
D = D(~ismember({D.name}, {'.', '..'})); % dir returns '.' and '..' (usually in first slot)

D = D(isfolder({D.name})); %only look at subdirectories
counter = 1; % to keep track of total iterations and row assignment

sent = strings(243,1);
condition = strings(243,1);
stepmorph = strings(243,1);
sumpitch = zeros(243,1);
maxpitch = zeros(243,1);
sumloud = zeros(243,1);
maxloud = zeros(243,1);

for ii = 1:numel(D) % iterate over every sound file
      currD = D(ii).name; % get the current subdirectory name
      sent_name = extractAfter(D(ii).name,'sent');
      cd(currD) % enter the subdirectory
      conditionFolders = dir;
      conditionFolders = conditionFolders(~ismember({conditionFolders.name}, {'.', '..'})); % dir returns '.' and '..' (usually in first slot)
      conditionFolders = conditionFolders(isfolder({conditionFolders.name}));
      for jj = 1:3
          currD = conditionFolders(jj).name;
          cd(currD)
          condition_name = extractAfter(conditionFolders(jj).name,'_');
          soundfiles = dir('*.wav'); % get all wav files in current subdirectory
%% Locate the neutral base
        for pp = 1:4
            current_soundname = extractAfter(soundfiles(pp).name,'00');
            if current_soundname == '1.wav'
                base_neutral = soundfiles(pp);
                [y1,fs] = audioread(base_neutral.name);
                break % exit loop once the neutral base is found
            end
        end
%% Analysis of Morphs 
        for kk = 1:4
            if kk == pp % if it is the neutral base, do not evaluate
                continue
            else
                [y2,~] = audioread(soundfiles(kk).name);
                [delta_pitch,delta_loud,~,~,~,~,~,~] = pitch_loud_diff(y1,y2,fs,params,false);
                
                current_morph_step = extractBetween(soundfiles(kk).name,'00','.wav');
                sent(counter,1) = sent_name;
                condition(counter,1) = condition_name;
                stepmorph(counter,1) = current_morph_step{1,1};
                sumpitch(counter,1) = nansum(delta_pitch);
                maxpitch(counter,1) = nanmax(delta_pitch);
                sumloud(counter,1) = nansum(delta_loud);
                maxloud(counter,1) = nanmax(delta_loud);
                
                counter = counter + 1;
            end  
        end
        disp(ii)
        cd('..') % return to submain directory
      end
      cd('..') % return to main directory
end
T = table(sent,condition,stepmorph,sumpitch,maxpitch,sumloud,maxloud);
writetable(T,'analysis.csv');