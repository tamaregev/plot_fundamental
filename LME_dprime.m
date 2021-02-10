%% Linear mixed effects analysis
pilotFolder = '/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/Results/Pilot/';
T=readtable([pilotFolder 'dprimes.xls']);

T.condition = categorical(T.condition);
C = categories(T.condition);
C2 = C([2,1,3]);
T.condition = reordercats(T.condition,C2);

formula = 'dprime ~ condition + (1+condition|subj)';
lme = fitlme(T,formula);

formula2 = 'dprime ~ condition + (1|subj)';
lme2 = fitlme(T,formula2);

compare(lme2,lme)

formula3 = 'dprime ~ condition';
lme3 = fitlme(T,formula3);

compare(lme3,lme)
