%% Constants
clear

% Length
m = 10^0; mm = 10^-3*m; um = 10^-6*m; nm = 10^-9*m;
% Time
s = 10^0; ps = 10^-12*s; fs = 10^-15*s;
% Energy
J = 10^0; uJ = 10^-6 * J; nJ = 10^-12 * J;
% Physics Constants
c = 299792458; eps0 = 8.854187817*10^-12;

% Run Parameters
songStr = 'Intro.mp3';

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

deltaKL = matlabFunction(...
    phi0o(lam1,l) + phi0o(lam2,l) - phi0e(lam3,l,theta) + (phi1*l) + (phi2*l),...
    'Vars',{lam1,lam2,lam3,theta,phi1,phi2,l},'File','bboShift');


nO = matlabFunction(nO,'Vars',{lam});
nE_Theta = matlabFunction(nE_Theta,'Vars',{lam,theta});

clear l theta lam lam1 lam2 lam3 phi1 phi2 phi0o phi0e 



%% Problem Definition

lams(1) = 1550 * nm;
lams(2) = lams(1);
lams(3) = 1/((1/lams(1))+(1/lams(2)));

omegas = 2*pi*c./lams;

crysLen = 3*mm;
thetaTST = [0 50];
[phis(1,:),songFs] = loadSongPhase(songStr,20,27,-3*pi/4,0,0);

phis(end,1) = 0;

figure(1);
plot(phis(1,:));
drawnow

N = length(phis(1,:));
phis(2,:) = zeros(N,1);

% Find crystal angle for phase matching
thetaZero = findTheta(lams,thetaTST,[-3*pi/4 0],crysLen);

pEnergy = 1 * uJ / 2;
pDur = 250 * fs;
pRad = 100 * um;

pInten = ((2*sqrt(log(16)))/((pi^(3/2))*(pRad^2)*(pDur)))*pEnergy;
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
        tRange, conds);
    ampVals(ii) = Y(end,3);
    timeVals(ii) = toc;
    if mod(ii,round(N/100)) == 0
        disp([num2str( round((ii/N)*100) ),'% complete']);
    end
end

disp(['Each iteration took: ', num2str(mean(timeVals)*1000,'%.2f'),' ms'])
disp(['The total time was: ', num2str(sum(timeVals)/60,'%.2f'),' min'])
ampVals = convAmp(ampVals);

figure(2);
plot(T,abs(Y).^2);
drawnow

figure(3);
plot(1:length(ampVals),ampVals);
drawnow

player = audioplayer(ampVals,songFs);
play(player)


function Sys = myODEs()
syms a1(t) a2(t) a3(t) w1 w2 w3 dNL c l1 l2 l3 th0 ph1 ph2 Y t
ode1 = diff(a1) == ((2i*w1*dNL)/c)*a3*conj(a2)*...
    exp(1i*bboShift(l1,l2,l3,th0,ph1,ph2,t));
ode2 = diff(a2) == ((2i*w2*dNL)/c)*a3*conj(a1)*...
    exp(1i*bboShift(l1,l2,l3,th0,ph1,ph2,t));
ode3 = diff(a3) == ((2i*w3*dNL)/c)*a1*a2*...
    exp(1i*bboShift(l1,l2,l3,th0,ph1,ph2,t));
[ODE,~] = odeToVectorField(ode1, ode2, ode3);
Sys = matlabFunction(ODE, 'Vars', {t, Y, w1, w2, w3, dNL, c, l1, l2, l3, th0, ph1, ph2});
end


function y = findTheta(lams,thetaTst,phis,crysLen)

y = fzero(@fun,thetaTst,optimset('Display','off'));

    function val = fun(theta)
        val = bboShift(lams(1),lams(2),lams(3),theta,phis(1),phis(2),crysLen);
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

function amp = convAmp(amp)

amp = abs(amp).^2;

amp = (( (amp - min(amp)) / (max(amp)-min(amp)) ) * 2) - 1;

end





