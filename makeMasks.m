function [mask_all, combinations] = makeMasks(folder,nMasks,nStim)
%MAKEMASKS makes masks that are the mean of soundfiles in the input folder.
% Both a mask that is the mean of all stimuli (mask_all) is calucalted and nMask
% that consists of random combinations of nStim out of the whole pool. 
% the masks are saved into a folder named Masks within the input folder.

if ~exist([folder filesep 'Masks'],'dir')
    mkdir(folder,'Masks')
end

Files = dir([folder filesep '*.wav']);
lengths=nan(length(Files),1);
clear y
for i=1:length(Files)
    [y{i},fs]=audioread([folder filesep Files(i).name]);
    lengths(i) = length(y{i});
end
minlen = min(lengths);

y=nan(minlen,length(Files));

for i=1:length(Files)
    [y1,fs]=audioread([folder filesep Files(i).name]);
    y(:,i)=y1(1:minlen);
end
mask_all = mean(y,2);
mask_all = hann_fade(mask_all,100,fs);

% soundsc(mask_all,fs)
% pause
audiowrite([folder filesep 'Masks' filesep 'mask_all.wav'],mask_all,fs);

% random groups of nStim
combinations = nan(nMasks, nStim);

for ig = 1:nMasks
    whichFiles = randsample(length(Files),nStim);
    combinations(ig,:) = whichFiles;
    
    lengths = nan(length(whichFiles),1);
    clear y
    i=1;
    for fi = whichFiles'
        [y{i},fs] = audioread([folder filesep Files(fi).name]);
        lengths(i) = length(y{i});
        i=i+1;
    end
    minlen = min(lengths);

    y=nan(minlen,length(whichFiles));
    i=1;
    for fi=whichFiles'
        [y1,fs]=audioread([folder filesep Files(fi).name]);
        y(:,i)=y1(1:minlen);
        i=i+1;
    end
    mask = mean(y,2);
    mask = hann_fade(mask,100,fs);

%     soundsc(mask,fs)
%     pause

    audiowrite([folder filesep 'Masks' filesep 'mask_' num2str(ig) '.wav'],mask,fs);
end
%
save([folder filesep 'Masks' filesep 'cominations'],'combinations');

end

