% catmodel - [Zilany, Bruce, Nelson and Carney] Auditory Nerve Model
%
%     vihc = catmodel_IHC(pin,CF,nrep,binwidth,reptime,cohc,cihc);
%     [synout,psth] = catmodel_Synapse(vihc,CF,nrep,binwidth,reptime,fibertype,implnt);
%
% vihc is the inner hair cell (IHC) potential (in volts)
% synout is the synapse output
% psth is the peri-stimulus time histogram 
%
% pin is the input sound wave in Pa sampled at the appropriate sampling rate (see instructions below)
% CF is the characteristic frequency of the fiber in Hz
% nrep is the number of repetitions for the psth
% binwidth is the binsize in seconds, i.e., the reciprocal of the sampling rate (see instructions below)
% reptime is the time between stimulus repetitions in seconds - NOTE should be equal to or longer than the duration of pin
% cohc is the OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
% cihc is the IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
% fibertype is the type of the fiber based on spontaneous rate in spikes/s - "1" for Low SR; "2" for Medium SR; "3" for High SR
% implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
% 
%
% Fore example,
%
%    vihc = catmodel_IHC(pin,1e3,10,1/100e3,0.200,1,1);
%    [synout,psth] = catmodel_Synapse(vihc,1e3,10,1/100e3,0.200,3,0);
%
% models a normal fiber of high spontaneous rate (normal OHC & IHC function) with a CF of 1 kHz, 
% for 10 repititions and a sampling rate of 100kHz, for a repetition duration of 200 ms, and
% with approximate implementation of the power-law functions in the synapse model.
%
%
% NOTE ON SAMPLING RATE:-
% Since version 2 of the code, it is possible to run the model at a range
% of sampling rates between 100 kHz and 500 kHz.
% It is recommended to run the model at 100 kHz for CFs up to 20 kHz, and
% at 200 kHz for CFs up to 40 kHz.