function analyze_allsubjects_data()
% 
% analyze_allsubjects_data() reads trial data from the current directory
% and graphs either average performance or individual performance for all
% conditions
%
% Note: comment out"Plot averages" or "Plot individual data",
% depending on what you want to graph OR use "Run Section"
%
% August 2020, Julie Meng, initial version
%% Read data files
addpath(genpath('Users/julie/Downloads/ProsodyUROP/CodedExperiment/simpledata'));
Files = dir('*.csv');

FileNames = string.empty;
for k=1:length(Files)
   FileNames = [FileNames Files(k).name];
end

for k=1:length(FileNames)
    disp(FileNames(k))
end


%% Create tables

Q_all = [0 0 0 0];
F_all = [0 0 0 0];
No_all = [0 0 0 0];
morphs = [1 2 3 4];

 for index = 1:length(FileNames)
    disp(strcat('simpledata/', FileNames(index)));
    T = readtable(strcat('simpledata/', FileNames(index)));

    %figure out level of morph
    T.morph = max(T.Level1, T.Level2);

    %select only critical
    T = sortrows(T, 'Critical_Filler');
    C = T(strcmp(T.Critical_Filler, 'critical'), :);

    %table of critical Q
    Q = C(strcmp(C.Condition, 'Q'), :);

    %make table of morph and percent differences
    Q_1 = height(Q((Q.morph == 1) & strcmp(Q.SOrD, 'D'),:))/height(Q((Q.morph == 1),:));
    Q_2 = height(Q((Q.morph == 2) & strcmp(Q.SOrD, 'D'),:))/height(Q((Q.morph == 2),:));
    Q_3 = height(Q((Q.morph == 3) & strcmp(Q.SOrD, 'D'),:))/height(Q((Q.morph == 3),:));
    Q_4 = height(Q((Q.morph == 4) & strcmp(Q.SOrD, 'D'),:))/height(Q((Q.morph == 4),:));
    
    Q_perc = [Q_1 Q_2 Q_3 Q_4];


    %table of critical F
    F = C(strcmp(C.Condition, 'F'), :);

    %make table of morph and percent differences
    F_1 = height(F((F.morph == 1) & strcmp(F.SOrD, 'D'),:))/height(F((F.morph == 1),:));
    F_2 = height(F((F.morph == 2) & strcmp(F.SOrD, 'D'),:))/height(F((F.morph == 2),:));
    F_3 = height(F((F.morph == 3) & strcmp(F.SOrD, 'D'),:))/height(F((F.morph == 3),:));
    F_4 = height(F((F.morph == 4) & strcmp(F.SOrD, 'D'),:))/height(F((F.morph == 4),:));
    F_perc = [F_1 F_2 F_3 F_4];

    %table of critical No
    No = C(strcmp(C.Condition, 'No'), :);

    %make table of morph and percent differences
    No_1 = height(No((No.morph == 1) & strcmp(No.SOrD, 'D'),:))/height(No((No.morph == 1),:));
    No_2 = height(No((No.morph == 2) & strcmp(No.SOrD, 'D'),:))/height(No((No.morph == 2),:));
    No_3 = height(No((No.morph == 3) & strcmp(No.SOrD, 'D'),:))/height(No((No.morph == 3),:));
    No_4 = height(No((No.morph == 4) & strcmp(No.SOrD, 'D'),:))/height(No((No.morph == 4),:));
    No_perc = [No_1 No_2 No_3 No_4];
        
     Q_all = cat(1, Q_all, Q_perc);
     F_all = cat(1, F_all, F_perc);
     No_all = cat(1, No_all, No_perc);
    
 end
 
 Q_all(1,:) = [];
%  disp(Q_all);
 
 F_all(1,:) = [];
%  disp(F_all);
 
 No_all(1,:) = [];
%  disp(No_all);
    
    
%% Plot averages

hold on
 
%   plot(morphs, mean(Q_all), '-o');
%   plot(morphs, mean(F_all), '-o');
%   plot(morphs, mean(No_all), '-o');
%   
  errorbar(morphs, mean(Q_all), std(Q_all)/sqrt(length(Q_all)), '-o', 'LineWidth', 1.5);
  errorbar(morphs, mean(F_all), std(F_all)/sqrt(length(F_all)), '-o', 'LineWidth', 1.5);
  errorbar(morphs, mean(No_all), std(No_all)/sqrt(length(No_all)), '-o', 'LineWidth', 1.5);
  
  t = title('Average performance across all subjects');
  t.FontSize = 18;
  ylab = ylabel('% Different');
  ylab.FontSize = 18;
  xlab = xlabel('Morph #');
  xlab.FontSize = 18;
  xticks([1 2 3 4])
  yticks([0 .2 .4 .6 .8 1])
  yticklabels({'0' '20' '40' '60' '80' '100'})
  
 hold off
 
 lgd = legend('Q', 'F', 'No', 'Location', 'southeast');
 lgd.FontSize = 15;
 lgd.Title.String = 'Condition';
 
 %% Plot individual data
 
%  plot(morphs, Q_all, '-o');

 subplot(1,3,1);
 plot(morphs, Q_all, '-o', 'Color', [0    0.4470    0.7410], 'LineWidth', 1.5);
 q_t = title('Individuals - Q');
 q_t.FontSize = 16;
 xlabel('Morph #', 'FontSize', 14)
 ylabel('% Different', 'FontSize', 14)
 xticks([1 2 3 4])
 yticks([0 .2 .4 .6 .8 1])
 yticklabels({'0' '20' '40' '60' '80' '100'})
 
 subplot(1,3,2);
 plot(morphs, F_all, '-o', 'Color', [0.8500    0.3250    0.0980], 'LineWidth', 1.5)
 f_t = title('Individuals - F');
 f_t.FontSize = 16;
 xlabel('Morph #', 'FontSize', 14)
 ylabel('% Different', 'FontSize', 14)
 xticks([1 2 3 4])
 yticks([0 .2 .4 .6 .8 1])
 yticklabels({'0' '20' '40' '60' '80' '100'})

 
 subplot(1,3,3);
 plot(morphs, No_all, '-o', 'Color', [0.9290    0.6940    0.1250], 'LineWidth', 1.5)
 no_t = title('Individuals - No');
 no_t.FontSize = 16;
 xlabel('Morph #', 'FontSize', 14)
 ylabel('% Different', 'FontSize', 14)  
 xticks([1 2 3 4])
 yticks([0 .2 .4 .6 .8 1])
 yticklabels({'0' '20' '40' '60' '80' '100'})
     
 end
