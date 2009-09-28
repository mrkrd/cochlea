function [Cohc,Cihc,OHC_Loss]=fitaudiogram(FREQUENCIES,dBLoss,Dsd_OHC_Loss)

% FITAUDIOGRAM Gives values of Cohc and Cihc that produce a desired
% threshold shift for the cat auditory-periphery model of Zilany and Bruce
% (J. Acoust. Soc. Am. 2007).
%
% [Cohc,Cihc,OHC_Loss]=fitaudiogram(FREQUENCIES,dBLoss,Dsd_OHC_Loss)\
%
% The output variables are arrays with values corresponding to each
% frequency in the input array FREQUENCIES.
%
% Cohc is the outer hair cell (OHC) impairment factor; a value of 1
%      corresponds to normal OHC function and a value of 0 corresponds to
%      total impairment of the OHCs.
%
% Cihc is the inner hair cell (IHC) impairment factor; a value of 1
%      corresponds to normal IHC function and a value of 0 corresponds to
%      total impairment of the IHC.
%
% OHC_Loss is the threshold shift in dB that is attributed to OHC
%      impairment (the remainder is produced by IHC impairment).
%
%
% FREQUENCIES is an array of frequencies in Hz for which you have audiogram data.
%
% dBLoss is an array of threshold shifts in dB (for each frequency in FREQUENCIES).
%
% Dsd_OHC_Loss is an optional array giving the desired threshold shift in
%      dB that is caused by the OHC impairment alone (for each frequency
%      in FREQUENCIES).
%      If this array is not given, then the default desired threshold shift
%      due to OHC impairment is 2/3 of the entire threshold shift at each
%      frequency.  This default is consistent with the effects of acoustic
%      trauma in cats (see Bruce et al., JASA 2003, and Zilany and Bruce,
%      JASA 2007) and estimated OHC impairment in humans (see Moore,
%      Glasberg & Vickers, JASA 1999).
%
% © M. S. A. Zilany and I. C. Bruce (ibruce@ieee.org), March - April 2007


load THRESHOLD_TUNING_ALL;
% Variables are
% CF: 125 Hz to 10 kHz                   [1*37]
% CIHC: varies from 1.0 to 0.0001        [1*56]
% COHC: varies from 1.0 to 0             [1*56]
% THR : absolute thresholds              [37*56*56]

for k = 1:length(THR(:,1,1))
    dBShift(k,:,:)= THR(k,:,:) - THR(k,1,1);
end

if nargin<3, Dsd_OHC_Loss = 2/3*dBLoss;
end;

for m = 1:length(FREQUENCIES)
    [W,N] = min(abs(CF-FREQUENCIES(m))); n = N(1);
    
    if Dsd_OHC_Loss(m)>dBShift(n,1,end)
        Cohc(m) = 0;
    else
        [a,idx]=sort(abs(squeeze(dBShift(n,1,:))-Dsd_OHC_Loss(m)));
        Cohc(m)=COHC(idx(1));
    end
    OHC_Loss(m) = interp1(COHC,squeeze(dBShift(n,1,:)),Cohc(m),'nearest');
    [mag,ind] = sort(abs(COHC-Cohc(m)));
    
    Loss_IHC(m) = dBLoss(m)-OHC_Loss(m);
    
    if dBLoss(m)>dBShift(n,end,ind(1))
        Cihc(m) = 0;
    else
        [c,indx]=sort(abs(squeeze(dBShift(n,:,ind(1)))-dBLoss(m)));
        Cihc(m)=CIHC(indx(1));
    end
    
end
    
    