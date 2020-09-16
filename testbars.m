% x_sum = categorical(tableb.Properties.VariableNames(3:2:5));
% x_sum = reordercats(x_sum,tableb.Properties.VariableNames(3:2:5));
% 
% x_max = categorical(tableb.Properties.VariableNames(4:2:6));
% x_max = reordercats(x_max,tableb.Properties.VariableNames(4:2:6));
stde_den = sqrt(18);
% stde_den2 = sqrt(9);
%% Setting variables - critical sentences
mp_no1to2 = std(T.maxpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
mp_f1to2 = std(T.maxpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
mp_q1to2 = std(T.maxpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
y_mp1to2 = [mp_q1to2,mp_f1to2,mp_no1to2];
y_mp_stderror1to2 = y_mp1to2./stde_den;

mp_no1to3 = std(T.maxpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
mp_f1to3 = std(T.maxpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
mp_q1to3 = std(T.maxpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
y_mp1to3 = [mp_q1to3,mp_f1to3,mp_no1to3];
y_mp_stderror1to3 = y_mp1to3./stde_den;

mp_no1to4 = std(T.maxpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
mp_f1to4 = std(T.maxpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
mp_q1to4 = std(T.maxpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
y_mp1to4 = [mp_q1to4,mp_f1to4,mp_no1to4];
y_mp_stderror1to4 = y_mp1to4./stde_den;

Ymp = [mp_q1to2,mp_q1to3,mp_q1to4;mp_f1to2,mp_f1to3,mp_f1to4;mp_no1to2,mp_no1to3,mp_no1to4];
Ymp_se = [y_mp_stderror1to2;y_mp_stderror1to3;y_mp_stderror1to4];
%%
sp_no1to2 = std(T.sumpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
sp_f1to2 = std(T.sumpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
sp_q1to2 = std(T.sumpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
y_sp1to2 = [sp_q1to2,sp_f1to2,sp_no1to2];
y_sp_stderror1to2 = y_sp1to2./stde_den;

sp_no1to3 = std(T.sumpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
sp_f1to3 = std(T.sumpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
sp_q1to3 = std(T.sumpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
y_sp1to3 = [sp_q1to3,sp_f1to3,sp_no1to3];
y_sp_stderror1to3 = y_sp1to3./stde_den;

sp_no1to4 = std(T.sumpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
sp_f1to4 = std(T.sumpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
sp_q1to4 = std(T.sumpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
y_sp1to4 = [sp_q1to4,sp_f1to4,sp_no1to4];
y_sp_stderror1to4 = y_sp1to4./stde_den;

Ysp = [sp_q1to2,sp_q1to3,sp_q1to4;sp_f1to2,sp_f1to3,sp_f1to4;sp_no1to2,sp_no1to3,sp_no1to4];
Ysp_se = [y_sp_stderror1to2;y_sp_stderror1to3;y_sp_stderror1to4];
%%
ml_no1to2 = std(T.maxloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
ml_f1to2 = std(T.maxloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
ml_q1to2 = std(T.maxloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
y_ml1to2 = [ml_q1to2,ml_f1to2,ml_no1to2];
y_ml_stderror1to2 = y_ml1to2./stde_den;

ml_no1to3 = std(T.maxloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
ml_f1to3 = std(T.maxloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
ml_q1to3 = std(T.maxloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
y_ml1to3 = [ml_q1to3,ml_f1to3,ml_no1to3];
y_ml_stderror1to3 = y_ml1to3./stde_den;

ml_no1to4 = std(T.maxloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
ml_f1to4 = std(T.maxloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
ml_q1to4 = std(T.maxloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
y_ml1to4 = [ml_q1to4,ml_f1to4,ml_no1to4];
y_ml_stderror1to4 = y_ml1to4./stde_den;

Yml = [ml_q1to2,ml_q1to3,ml_q1to4;ml_f1to2,ml_f1to3,ml_f1to4;ml_no1to2,ml_no1to3,ml_no1to4];
Yml_se = [y_ml_stderror1to2;y_ml_stderror1to3;y_ml_stderror1to4];
%%
sl_no1to2 = std(T.sumloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
sl_f1to2 = std(T.sumloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
sl_q1to2 = std(T.sumloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),2)));
y_sl1to2 = [sl_q1to2,sl_f1to2,sl_no1to2];
y_sl_stderror1to2 = y_sl1to2./stde_den;

sl_no1to3 = std(T.sumloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
sl_f1to3 = std(T.sumloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
sl_q1to3 = std(T.sumloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),3)));
y_sl1to3 = [sl_q1to3,sl_f1to3,sl_no1to3];
y_sl_stderror1to3 = y_sl1to3./stde_den;

sl_no1to4 = std(T.sumloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
sl_f1to4 = std(T.sumloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
sl_q1to4 = std(T.sumloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[1:18]) & ismember(str2double(T.stepmorph),4)));
y_sl1to4 = [sl_q1to4,sl_f1to4,sl_no1to4];
y_sl_stderror1to4 = y_sl1to4./stde_den;

Ysl = [sl_q1to2,sl_q1to3,sl_q1to4;sl_f1to2,sl_f1to3,sl_f1to4;sl_no1to2,sl_no1to3,sl_no1to4];
Ysl_se = [y_sl_stderror1to2;y_sl_stderror1to3;y_sl_stderror1to4];
%% Setting variables - filler sentences
% mp_no2 = T.maxpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[19:27]));
% mp_f2 = T.maxpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[19:27]));
% mp_q2 = T.maxpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[19:27]));
% y_mp2 = [std(mp_q2),std(mp_f2),std(mp_no2)];
% y_mp_stderror2 = y_mp2./stde_den2;
% 
% sp_no2 = T.sumpitch(strcmp(T.condition,'No') & ismember(str2double(T.sent),[19:27]));
% sp_f2 = T.sumpitch(strcmp(T.condition,'F') & ismember(str2double(T.sent),[19:27]));
% sp_q2 = T.sumpitch(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[19:27]));
% y_sp2 = [std(sp_q2),std(sp_f2),std(sp_no2)];
% y_sp_stderror2 = y_sp2./stde_den2;
% 
% ml_no2 = T.maxloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[19:27]));
% ml_f2 = T.maxloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[19:27]));
% ml_q2 = T.maxloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[19:27]));
% y_ml2 = [std(ml_q2),std(ml_f2),std(ml_no2)];
% y_ml_stderror2 = y_ml2./stde_den2;
% 
% sl_no2 = T.sumloud(strcmp(T.condition,'No') & ismember(str2double(T.sent),[19:27]));
% sl_f2 = T.sumloud(strcmp(T.condition,'F') & ismember(str2double(T.sent),[19:27]));
% sl_q2 = T.sumloud(strcmp(T.condition,'Q') & ismember(str2double(T.sent),[19:27]));
% y_sl2 = [std(sl_q2),std(sl_f2),std(sl_no2)];
% y_sl_stderror2 = y_sl2./stde_den2;
%%
figure

subplot(2,2,1)
bar(Ymp)
set(gca, 'XTickLabel', {'Q','F','No'});
hold on
ngroups = size(Ymp,1);
nbars = size(Ymp,2);
groupwidth = min(.8,nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/ (2*nbars);
    er = errorbar(x,Ymp(:,i),Ymp_se(:,i),'.');
    er.LineStyle = 'none';
    er.LineWidth = 2;
    er.Color = 'k';
end
hold off
title('Maximum Difference in Pitch')
ylabel('Hz')

subplot(2,2,2)
bar(Ysp)
set(gca, 'XTickLabel', {'Q','F','No'});
hold on
ngroups = size(Ysp,1);
nbars = size(Ysp,2);
groupwidth = min(.8,nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/ (2*nbars);
    er = errorbar(x,Ysp(:,i),Ysp_se(:,i),'.');
    er.LineStyle = 'none';
    er.LineWidth = 2;
    er.Color = 'k';
end
hold off
title('Sum of the Difference in Pitch')
ylabel('Hz')

subplot(2,2,3)
bar(Yml)
set(gca, 'XTickLabel', {'Q','F','No'});
hold on
ngroups = size(Yml,1);
nbars = size(Yml,2);
groupwidth = min(.8,nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/ (2*nbars);
    er = errorbar(x,Yml(:,i),Yml_se(:,i),'.');
    er.LineStyle = 'none';
    er.LineWidth = 2;
    er.Color = 'k';
end
hold off
title('Maximum Difference in Loudness')
ylabel('SPL (dB)')

subplot(2,2,4)
bar(Ysl)
set(gca, 'XTickLabel', {'Q','F','No'});
hold on
ngroups = size(Ysl,1);
nbars = size(Ysl,2);
groupwidth = min(.8,nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/ (2*nbars);
    er = errorbar(x,Ysl(:,i),Ysl_se(:,i),'.');
    er.LineStyle = 'none';
    er.LineWidth = 2;
    er.Color = 'k';
end
title('Sum of the Difference in Loudness')
ylabel('SPL (dB)')
legend('step1 to step2','step1 to step3','step1 to step4')
