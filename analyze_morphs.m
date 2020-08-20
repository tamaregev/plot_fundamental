% current folder should be main folder with subdirectories
params.th_f0score = 0.75; %threshold f0 power for plotting f0
params.th_df = 95; %Hz maximal f0 jump for plotting f0
params.conv = 5;
params.my_cal_factor = 1.0;  %the value for your system to convert the WAV into Pascals

D = dir; % A is a struct ... first elements are '.' and '..' used for navigation.
D = D(~ismember({D.name}, {'.', '..'})); % dir returns '.' and '..' (usually in first slot)
                                            % , filter them out for eval.

D = D(isfolder({D.name})); %only look at subdirectories
% varheaders = ["sent","condition","sum_pitch_diff","max_pitch_diff",...
          % "sum_loud_diff","max_loud_diff"]; % table headers
counter = 1; % to keep track of total iterations and row assignment

sent = strings(81,1);
condition = strings(81,1);
sumpitch = zeros(81,1);
maxpitch = zeros(81,1);
sumloud = zeros(81,1);
maxloud = zeros(81,1);

for ii = 1:numel(D) % iterate over every sound file
      currD = D(ii).name; % get the current subdirectory name
      sent_name = extractAfter(D(ii).name,'Sent');
      cd(currD) % enter the subdirectory
      soundfiles = dir('*.wav'); % get all wav files in current subdirectory
% %% Identify which wav file is the neutral base
%       pattern = 'N.wav';      
%       for jj = 1:numel(soundfiles)
%         currenttest = endsWith(soundfiles(jj).name, pattern);
%         if currenttest
%             neutralbase = soundfiles(jj); % store neutral base file for comparison
%             baseelementnumber = jj;
%         end
%       end
%% Locate the neutral base
%       base_neutral = soundfiles(2);
%       [y1,fs] = audioread(base_neutral.name);
      for pp = 1:4
          current_soundname = extractBetween(soundfiles(pp).name,'_','.wav');
          if current_soundname{1,1} == 'N'
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
              [delta_pitch,delta_loud, p1, p2, t_p, l1, l2, t_l] = pitch_loud_diff(y1,y2,fs,params,false);
              
              sent(counter,1) = sent_name;
              cellarray_cond = extractBetween(soundfiles(kk).name,'_','.wav');
              condition(counter,1) = cellarray_cond{1,1};
              sumpitch(counter,1) = nansum(delta_pitch);
              maxpitch(counter,1) = nanmax(delta_pitch);
              sumloud(counter,1) = nansum(delta_loud);
              maxloud(counter,1) = nanmax(delta_loud);
              
              counter = counter + 1;
          end  
      end
      disp(ii)
      cd('..') % return to main directory
end
T = table(sent,condition,sumpitch,maxpitch,sumloud,maxloud);
writetable(T,'analysis.csv');