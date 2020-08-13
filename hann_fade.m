function output = hann_fade(input,ramp_dur_ms,fs)
%HANN_FADE applies a (half) hanning window fade in and fade out (hanning provides a smoother function than a linear ramp) 

% Tamar Regev, Aug 2020
% Adapted from Malinda's Hanning Window Code (hann)


if nargin < 3
   fs = 41000;
end
   
output = input;

stim_dur_smp = length(input);
ramp_dur_smp = floor(ramp_dur_ms * fs / 1000);
if (stim_dur_smp < 2*ramp_dur_smp)
 error('Ramps cannot be longer than the stimulus duration')
end


win = hanning(ramp_dur_smp*2);

%fade in (ramp up) initial part:

output(1:ramp_dur_smp) = win(1:ramp_dur_smp) .* output(1:ramp_dur_smp);

%fade out (ramp down) last part:
output(end-ramp_dur_smp+1:end) = win(ramp_dur_smp+1:ramp_dur_smp*2) .* output(end-ramp_dur_smp+1:end);
