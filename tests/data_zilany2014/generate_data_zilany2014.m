clear all;


%%% Model parameters
cf = 1e3;
cohc = 1.0;
cihc = 1.0;
species = 1; % 1: cat; 2: human with Shera et al. tuning; 3: for human with Glasberg & Moore tuning
noise_type = 0;
fiber_type = 3; % 3: HSR
implnt = 0; % 0: approximate; 1: actual
nrep = 1;


%%% Stimulus
fs = 100e3;
tmax = 50e-3;

sound = zeros(1, fs*tmax);
sound(10e-3*fs) = 1;


%%% Run simmulation
vihc = model_IHC(sound, cf, nrep, 1/fs, tmax, cohc, cihc, species);
[meanrate,varrate,psth] = model_Synapse(vihc, cf, nrep, 1/fs, fiber_type, noise_type, implnt);


%%% Save results
save('data_zilany2014.mat', 'fs', 'cf', 'sound', 'vihc', 'meanrate')


%%% Plots
subplot(4,1,1)
plot(sound)

subplot(4,1,2)
plot(vihc)

subplot(4,1,3)
plot(meanrate)

subplot(4,1,4)
plot(psth)
