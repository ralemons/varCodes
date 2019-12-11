%% Constants
clear

% Length
m = 10^0; mm = 10^-3*m; um = 10^-6*m; nm = 10^-9*m;
% Time
s = 10^0; ps = 10^-12*s; fs = 10^-15*s; Hz = 1/s;
% Energy
J = 10^0; uJ = 10^-6 * J; nJ = 10^-9 * J;
% Physics Constants
c = 299792458; eps0 = 8.854187817*10^-12;

% Run Parameters
% songStr = 'Intro.mp3';
fileStr = '~/Documents/GitHub/varCodes/nonlinODEsolve/bboShift1.m';

%% Equations

syms l theta lam lam1 lam2 lam3 phi1 phi2

nO(lam) = sqrt( 2.7359 + 0.01878/((lam/10^-6)^2 - 0.01822)...
    - 0.01354 * (lam/10^-6)^2 );
nE(lam) = sqrt( 2.3753 + 0.01224/((lam/10^-6)^2 - 0.01667)...
    - 0.01516 * (lam/10^-6)^2 );
dNL = 2.01 * 10^-12;

nE_Theta(lam,theta) = sqrt( 1/ (cosd(theta)^2/nO(lam)^2 + sind(theta)^2/nE(lam)^2) );

phi0o(lam,l) = (2*pi/lam) * nO(lam) * l;
phi0e(lam,l,theta) = (2*pi/lam) * nE_Theta(lam,theta) * l;

if ~isfile(fileStr)
deltaKL = matlabFunction(...
    phi0o(lam1,l) + phi0o(lam2,l) - phi0e(lam3,l,theta) + (phi2) + (phi2),...
    'Vars',{lam1,lam2,lam3,theta,phi1,phi2,l},'File',fileStr);
end

nO = matlabFunction(nO,'Vars',{lam});
nE_Theta = matlabFunction(nE_Theta,'Vars',{lam,theta});

clear l theta lam lam1 lam2 lam3 phi1 phi2 phi0o phi0e 



%% Problem Definition

lams(1) = 1550 * nm;
lams(2) = lams(1);
lams(3) = 1/((1/lams(1))+(1/lams(2)));

omegas = 2*pi*c./lams;

crysLen = 2*mm;
thetaTST = [0 50];

useSin = 0;

if exist('songStr','var')
    
    [phis(1,:),sigFs] = loadSongPhase(songStr,24,25,-9*pi/4,0,0);
    phis(1,end) = 0;
    
elseif useSin
    
    sigFs = 100 * Hz;
    
    len = 20;
    scale = 3*pi;
    shift = 0;
    x = 0:1/sigFs:len;
    
    phis(1,:) = (scale*sin(2*pi*x/len))-shift;
    phis(1,end) = 0;
    
else
    
    Str = 'This is my test'; %#ok<*UNRCH>
    Bin = reshape(dec2bin(Str,8).',1,[])-'0';
    Bin = interpBin(Bin,2000);
    
    sigFs = 100 * Hz;
    scale = pi/2;
    shift = 0;
    
    phis(1,:) = (scale*Bin)+shift;
    
end

% phis(1,:) = 0;

% figure(1);
% plot(phis(1,:));
% drawnow

N = length(phis(1,:));


phis(2,:) = zeros(N,1);

clear saveData
M = 1;
saveData = struct('indVar','matching angle','var',zeros(1,M),'sigAmp',zeros(N,M),'T',{cell(1,M)},'Y',{cell(1,M)});

options = odeset('RelTol',1e-4,'AbsTol',1e-6);

% Find crystal angle for phase matching
for jj = 1:M
thetaZero = findTheta(lams,thetaTST,[-scale/2 0],crysLen);
saveData.var(1,jj) = thetaZero;
    
pEnergy = 1 * uJ;
% saveData.var(1,jj) = pEnergy;

pDur = 250 * fs;
pRad = 100 * um;

pInten = ((2*sqrt(log(16)))./((pi^(3/2)).*(pRad.^2).*(pDur))).*pEnergy;
pField = sqrt( (pInten)./(2*c*eps0*nO(lams(1))) );

tRange = [0 crysLen];
conds = [pField; pField; 0];

odes = myODEs();

% N = 100;
timeVals = zeros(N,1);
ampVals = zeros(N,1);

for ii = 1:N
    tic
    [T,Y] = ode45(@(t,Y)odes(t,Y,omegas(1),omegas(2),omegas(3),...
        dNL,c,lams(1),lams(2),lams(3),thetaZero,phis(1,ii),phis(2,ii)),...
        tRange, conds,options);
    ampVals(ii) = Y(end,3);
    timeVals(ii) = toc;
    if mod(ii,round(N/100)) == 0
        disp([num2str( round((ii/N)*100) ),'% complete']);
    end
end

disp(['Each iteration took: ', num2str(mean(timeVals)*1000,'%.2f'),' ms'])
disp(['The total time was: ', num2str(sum(timeVals)/60,'%.2f'),' min'])

normAmp = 0;
ampVals = convAmp(ampVals,normAmp);

saveData.sigAmp(:,jj) = ampVals;
saveData.T{jj} = T;
saveData.Y{jj} = abs(Y).^2;

figure(2);
plot(T,abs(Y).^2);
drawnow

figure(3);
plot(1:length(ampVals),ampVals);
drawnow

if normAmp
    player = audioplayer(ampVals,sigFs); %#ok<TNMLP>
    play(player)
end

end

function Sys = myODEs()
syms a1(t) a2(t) a3(t) w1 w2 w3 dNL c l1 l2 l3 th0 ph1 ph2 Y t
ode1 = diff(a1) == ((2i*w1*dNL)/c)*a3*conj(a2)*...
    exp(1i*bboShift(l1,l2,l3,th0,ph1,ph2,t,1));
ode2 = diff(a2) == ((2i*w2*dNL)/c)*a3*conj(a1)*...
    exp(1i*bboShift(l1,l2,l3,th0,ph1,ph2,t,2));
ode3 = diff(a3) == ((2i*w3*dNL)/c)*a1*a2*...
    exp(1i*bboShift(l1,l2,l3,th0,ph1,ph2,t,3));
[ODE,~] = odeToVectorField(ode1, ode2, ode3);
Sys = matlabFunction(ODE, 'Vars', {t, Y, w1, w2, w3, dNL, c, l1, l2, l3, th0, ph1, ph2});
end


function y = findTheta(lams,thetaTst,phis,crysLen)

y = fzero(@fun,thetaTst,optimset('Display','off'));

    function val = fun(theta)
        val = bboShift(lams(1),lams(2),lams(3),theta,phis(1),phis(2),crysLen,3);
    end

end

function [phi,Fs] = loadSongPhase(str,startVal,endVal,scale,shift,plotFlag)

[y,Fs] = audioread(str);

startVal = startVal*Fs;
endVal = endVal*Fs;
if startVal == 0
    startVal = 1;
end


phi = y(startVal:endVal,1);
phi = (scale*((phi )./max(phi )))+shift;

if exist('plotFlag','var') && plotFlag == 1
    plot(linspace(startVal/Fs,endVal/Fs,length(phi)),phi)
end

end

function amp = convAmp(amp,normFlag)

amp = abs(amp).^2;

if normFlag
    if length(amp) > 1
        amp = (( (amp - min(amp)) / (max(amp)-min(amp)) ) * 2) - 1;
    else
        amp = amp/max(amp);
    end
end



end

function out = interpBin(bin,newLen)

N = length(bin);
multFac = ceil(newLen/N);
out = zeros(1,N*multFac);

for ii = 1:N
    out( 1, 1+(multFac*(ii-1)):ii*multFac) = bin(ii);
end

end



