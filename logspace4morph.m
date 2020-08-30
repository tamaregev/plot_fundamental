function [frequencies,morphs] = logspace4morph(f1,f2,Nf,Nm)
%LOGSPACE4MORPH calculates the equal logarithmic distance frequencies and
%morph steps from the two edge frequencies.

% INPUTS:
%   f1,f2 are the frequencies of the first and last morph, in Hz. E.g., f1=180;f2=360;
%   Nf is the desired number of stimuli, e.g., Nf=5;
%   Nm is the number of morphs you span in a linear space. E.g., Nm=100
%
% OUTPUTS:
%   frequencies are the desired Nf frequencies with equal logarithmic distances. 
%   morphs are the morph numbers, out of the Nm linearly spaced morphs,
%   that should be selected to obtain the frequencies. You should round
%   them to select morphs that will approximate the desired frequencies.

    c=nthroot(f2/f1,Nf-1);%this is the multiplicative factor between the frequencies: frequencies = [f1,c*f1,c^2*f1 ...]
    
    frequencies = nan(1,Nf);
    frequencies(1)=f1;
    for i=2:Nf
        frequencies(i) = c*frequencies(i-1);
    end
    x_f = @(f) 1+(f-f1)*(Nm-1)./(f2-f1); 
    morphs = x_f(frequencies);
    
end

