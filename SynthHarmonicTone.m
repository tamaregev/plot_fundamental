function [ data, fulltime ] = SynthHarmonicTone( f0, harmonics, phases, dur, ramp, silence, energy, fs )
%Tamar 10.18.2015
% synth tone with harmonic (all at equal amplitude), 
% determining the fundamental frequency and which harmonics to include

%Inputs - 
% f0 - fundamental frequency (not necessarily included)
% harmonics - which harmoonic numbers to include. E.g. [2:6, 10]
% phases - of each harmonic. should be same size as harmonics.
% dur, ramp, silence, fs - as in SynthPureTone
% energy - total energy of the sound. Will determine normalisation.
%
%
%
%   Important! notice that the total length of the sound will be
%                   dur + sum(silence)/1000
%
% OUTPUT:
%   data - a vector of sound. play using e.g. sound(data,fs)
%   fulltime - a time vector, e.g. as x axis for plotting the sound wave -
%   plot(fulltime,data)

    T = 1/fs;       
    t = 0 : T : dur; 
    nHarmonics = length(harmonics);
    
    if length(phases) == 1
        phases = linspace(phases,phases,nHarmonics);
    end
    
    F = harmonics * f0;
    A = sqrt(energy/nHarmonics);

    Y = zeros(length(t),nHarmonics);
    ts = repmat(t',1,nHarmonics);
    Fs = repmat(F,length(t),1);
    phasess = repmat(phases, length(t),1);
    Y = A * sin( 2 * pi .* ts .* Fs + phasess );
    data = squeeze(sum(Y,2)./nHarmonics)';

    %rampup
    N = ceil(ramp(1)*fs/1000);
    rampup = linspace(0, 1, N); 
    data(1:N)=data(1:N).*(rampup);
    %rampdown
    N = ceil(ramp(2)*fs/1000);
    rampdown = linspace(1, 0, N); 
    data(end-N+1:end)=data(end-N+1:end).*(rampdown);

    %add silence
    before = linspace(0,0,silence(1)*fs/1000);
    after = linspace(0,0,silence(2)*fs/1000);
    data = [before, data, after];
    fulltime = 0:T:(dur + sum(silence)/1000);

end

