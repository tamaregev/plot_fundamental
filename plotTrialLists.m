D = dir(['/Users/tamaregev/Dropbox/postdoc/Fedorenko/Prosody/Prosody-meaning/trial_lists/input_text_files']);
D = D(~ismember({D.name}, {'.', '..','.DS_Store'})); % dir returns '.' and '..' (usually in first slot)
Tl=cell(numel(D),1);
Sent = nan(99,numel(D));
Cond = cell(99,numel(D));
Crit_Fill = cell(99,numel(D));
Morph = nan(99,numel(D));
DMorph = nan(99,numel(D));

for il=1:numel(D)
    Tl{il}=readtable([D(il).folder filesep D(il).name]);
    Tl{il}.Properties.VariableNames{1} = 'sent';
    Tl{il}.Properties.VariableNames{2} = 'cond';
    Tl{il}.Properties.VariableNames{3} = 'morph';
    C_F = cell(height(Tl{il}),1);
    morph = nan(height(Tl{il}),1);
    d_morphs = nan(height(Tl{il}),1);
    for it=1:height(Tl{il})
        morph(it)=str2double(Tl{il}.morph{it}(3));
        sent = Tl{il}.sent(it);
        if sent <= 18
            C_F{it} = 'C';
            d_morphs(it) = morph(it)-1;
        else
            C_F{it} = 'F';
            d_morphs(it) = 0;
        end
    end
    Tl{il} = [Tl{il}(:,1:2), table(C_F), table(morph), table(d_morphs)];
    %for plotting:
    Sent(:,il) = Tl{il}.sent;
    Cond(:,il) = Tl{il}.cond;
    Crit_Fill(:,il) = Tl{il}.C_F;
    Morph(:,il) = Tl{il}.morph;
    DMorph(:,il) = Tl{il}.d_morphs;
end
Cond_coded = zeros(size(Cond));
Cond_coded(strcmp(Cond,'F'))=1;
Cond_coded(strcmp(Cond,'No'))=2;

Sent_th = zeros(size(Sent));
Sent_th(Sent<19)=1;

imagesc(Sent');title('sent')
imagesc(Sent_th');title('Crit or Fill')
imagesc(Cond_coded');title('condition')
imagesc(DMorph');title('morph step diff')

