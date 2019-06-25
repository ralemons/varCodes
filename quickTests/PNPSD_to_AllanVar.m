%% Code to generate Allan Varience Plots from Phase Noise PSD data
% I may integrate this with the dataPlot_PhaseNoise.m file later



%% Run this section to load the data
clear;
close all;

% Read in files
[dataIN,headers,fileStr] = dataImport_PhaseNoise(); %#ok<*SAGROW>
numFilesIn = size(dataIN,3);

% Turn off warning about number of peaks. QOL change, can be removed.
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Run this section to process the data and find peaks

dataOUT = dataIN;

optFreq = 1.934e14; % Optical Frequency in Hz


for ii = 1:numFilesIn
    
    loopBW(ii) = str2double(headers{88,2,ii});
    carrier.freq(ii) = str2double(headers{72,2,ii});
    carrier.volt(ii) = str2double(headers{74,2,ii});
    phaseMod(ii) = str2double(headers{78,2,ii});
    rmsJitter(ii) = str2double(headers{80,2,ii});
    intPhaseNoise.dBc(ii) = str2double(headers{81,2,ii});
    intPhaseNoise.rad(ii) =...
        sqrt( 2 * trapz(dataOUT(:,1,ii) , (10 .^ (dataOUT(:,2,ii) ./ 10) ) ) );
    jitter.rf(ii) = (1/(2*pi*carrier.freq(1))) * intPhaseNoise.rad(ii);
    jitter.opt(ii) = (1/(2*pi*optFreq)) * intPhaseNoise.rad(ii);
    
    minVal.dBc(ii) = min(dataOUT(:,2,ii));
    maxVal.dBc(ii) = max(dataOUT(:,2,ii));
    
    
end


%%%%% Create the rad^2/Hz data %%%%%
for ii = 1:numFilesIn
    
    dataOUT(:,3,ii) = 2 * (10 .^ (dataOUT(:,2,ii) ./ 10) );
    minVal.rad(ii) = min(dataOUT(:,3,ii));
    maxVal.rad(ii) = max(dataOUT(:,3,ii));
    
end

%% Setup and Plot Allan Varience

clear Data n tau tau_0 f_carrier S_y sigma modsigma modADEV ADEV

psd(:,1) = dataOUT(:,1);
psd(:,2) = dataOUT(:,3);

f_carrier = optFreq;

% n = 2.^(1:21);

tmp = logspace(-1,7,50);
tmp(tmp > 3*10^6) = [];
tau = 1./tmp;

% tau_0 = .001;
% tau = tau_0*n;
% tau = linspace(1,n(end),1000);
% tau = [0.005;0.01;0.02;0.04;0.08;0.16;0.32;0.64;1.28;2.56;5.12;10.24;20.48;40.96;81.92;163.84];
% n = [1;2;4;8;16;32;64;128;256;512;1024;2048;4096;8192;16384];

ADEV = zeros(length(tau),1);
modADEV = zeros(length(tau),1);


for jj=1:length(tau)

        S_y=psd(:,2); 
        
        sigma=S_y.*sin(pi*tau(jj)*psd(:,1)).^4./(pi*tau(jj)*psd(:,1)).^2;     %fluctuations (ADEV shows frequency fluctuations)
        ADEV(jj)=sqrt(2*trapz(psd(:,1),sigma));                 %Formula see Riehle page 58 and following
        
        %from Bernier_Theoretical analysis ... simplified formula:
        modsigma=S_y.*2.*sin(pi*tau(jj)*psd(:,1)).^6./(pi*tau(jj)*psd(:,1)).^4;
        modADEV(jj)=sqrt(trapz(psd(:,1),modsigma));
        
end


loglog(tau,modADEV,'-o')
xlabel('Frequency (Hz)','Interpreter','tex')
ylabel('Allan Deviation (Hz^{-1})','Interpreter','tex')
title('Allan Deviation of the 24 mrad PSD')
xlim([tau(1) tau(end)])

