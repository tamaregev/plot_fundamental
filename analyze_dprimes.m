function [T, Q_dprimes, F_dprimes, No_dprimes] = analyze_dprimes()
% will try to incorporate into other script later
%% read dprime files
pilotFolder = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/Results/Pilot/';
dFolder = [pilotFolder 'dprimes/'];
addpath(genpath(dFolder)); %file names are id#_dprimes.csv
Files = dir([dFolder '*.csv']);

FileNames = string.empty;
for k=1:length(Files)
   FileNames = [FileNames Files(k).name];
end

for k=1:length(FileNames)
    disp(FileNames(k))
end

%% calculate mean and se

Q_dprimes = [];
F_dprimes = [];
No_dprimes = [];
num_subjects = length(FileNames);

Ta = table;

for index = 1:length(FileNames)
    T = readtable(FileNames(index));
    Q_dprimes(end+1) = T{1,1};
    F_dprimes(end+1) = T{1,2};
    No_dprimes(end+1) = T{1,3};
end

mean_Q_dprimes = mean(Q_dprimes)
mean_F_dprimes = mean(F_dprimes)
mean_No_dprimes = mean(No_dprimes)

SE_Q_dprimes = std(Q_dprimes)/sqrt(num_subjects)
SE_F_dprimes = std(F_dprimes)/sqrt(num_subjects)
SE_No_dprimes = std(No_dprimes)/sqrt(num_subjects)



%% make d-prime table for stats
conditions = {'Q','F','No'};
Td = table;

for index = 1:length(FileNames)
    fn = FileNames{index};
    T = readtable(fn);
    for con = 1:length(conditions)
        Temp = table;
        Temp.subj = index;
        Temp.ID = fn(1:end-12);
        Temp.condition = conditions(con);
        Temp.dprime = T{1,con};    
        Td = [Td ; Temp];
    end
end
writetable(Td,[pilotFolder 'dprimes.xls'])
T=Td;
