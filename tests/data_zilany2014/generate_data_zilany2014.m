clear all;
CF    = 1.0e3; % CF in Hz;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
species = 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds

% PSTH parameters
nrep = 1;              % number of stimulus repetitions (e.g., 50);

sound = zeros(1, Fs*T);
sound(10e-3*Fs) = 1;

subplot(4,1,1)
plot(sound)

vihc = model_IHC(sound, CF, nrep, 1/Fs, T, cohc, cihc, species);


subplot(4,1,2)
plot(vihc)



[meanrate,varrate,psth] = model_Synapse(vihc, CF, nrep, 1/Fs, fiberType, implnt);

subplot(4,1,3)
plot(synout)

subplot(4,1,4)
plot(psth)
