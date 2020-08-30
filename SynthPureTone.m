function [ data, fulltime ] = SynthPureTone( f, dur, ramp, silence, reduction, fs, hanFlag )
%Tamar 11.3.2015
% synth pure sinusoidals (cos) with an amplitude ramp for prevention of clicks. 
%   Can also add silence at beginning and end for further prevetion of the clicks.  
%
% INPUTS:
%   
%   f - frequency in Hz
%   dur - duration in seconds
%   ramp - length of amplitude rise and decay at [beginning, end] in ms.
%   silence - length of silence added at [beginning, end] in ms.
%       (wavplay needs [0 15] in order to not click for some reason)
%   reduction - factor of reduction in amplitude
%   fs - sampling rate (usually 44,100)
%   
%   hanFlag - optional (deafault: false) if true, ramp will be applied by
%   hanning window instead of linear.
%
%   Important! notice that the total length of the sound will be
%                   dur + sum(silence)/1000
%
%   

% OUTPUT:
%   data - a vector of sound. play using e.g. sound(data,fs)
%   fulltime - a time vector, e.g. as x axis for plotting the sound wave -
%   plot(fulltime,data)

% Aug 2020 Tamar: added input HanFlag

if ~exist('hanFlag','var')
    hanFlag=false;
end
    
    T = 1/fs; %sampling period
    t = 0:T:dur;
    data = cos(2*pi*f*t);
    if hanFlag
        %rampup
        N = ceil(ramp(1)*fs/1000);
        win=hanning(N*2);
        rampup = win(1:N)'; 
        data(1:N)=data(1:N).*(rampup);
        %rampdown
        N = ceil(ramp(2)*fs/1000);
        win=hanning(N*2)';
        rampdown = win(end-N+1:end); 
        data(end-N+1:end)=data(end-N+1:end).*(rampdown);
        
    else
        %rampup
        N = ceil(ramp(1)*fs/1000);
        rampup = linspace(0, 1, N); 
        data(1:N)=data(1:N).*(rampup);
        %rampdown
        N = ceil(ramp(2)*fs/1000);
        rampdown = linspace(1, 0, N); 
        data(end-N+1:end)=data(end-N+1:end).*(rampdown);
    end
    %add silence
    before = linspace(0,0,silence(1)*fs/1000);
    after = linspace(0,0,silence(2)*fs/1000);
    data = [before, data, after];
    data=data*reduction;
    fulltime = 0:T:(dur + sum(silence)/1000);
end

