%% Modeling of dispersion and pulse stacking in birefringent crystals

clear;

%% Constants
s = 10^12; % Change this to set the units of time
ns = s*10^-9;
ps = s*10^-12;
fs = s*10^-15;

m = 10^3; % Change this to set the units of space
cm = m*10^-2;
mm = m*10^-3;
um = m*10^-6;
c = 3*10^8 * m/s;

J = 10^6;
uJ = J*10^-6;

crys = 'BBO'; % BBO or YVO4 or CAL
lambda = 0.257; % Enter lambda in um regardless of set space unit
figNum = 1; % number for the figure you want, use to not overide other plots
tau = 2*ps/1.665; % tau = FWHM/1.665
power = 5*uJ; % normalize the energy in the pulse to real numbers
dCrysBase = 3 * mm; % base length of crystals
numCrys = 4;

% sympref('FourierParameters',[(1/sqrt(2*pi)) 1]); % you can uncomment this to use the same fourier parameters as Mathematica

[ne,no] = nonLinCrysChoice(crys,lambda); % get refractive indicies
dn = ne - no; % difference in refractive index between crystal axis


%% Crystal Definition
thetaCrys = ones(numCrys,1) .* deg2rad(45); % angle of the crystals

thetaPols = zeros(numCrys,1);

phiCrys = ones(numCrys,1) .* 0; % precise tunning of crystals
dCrys = ones(numCrys,1) .* dCrysBase .* 2.^(0:(numCrys-1))'; % generate larger crystals


%% Function Definitions

% Symbolic functions in matlab are kinda odd. I think this could be
% improved if I used anonymous functions but that is for a different time.
% The general idea here is that I create symbolic functions in the
% definition and all of the variables in them are also made. The *Han
% varibles are functions like matlab normally uses that you can plug in
% usual grids of points to get respective outputs.


syms f(t) F(w) R(theta) Pol(theta) BP(phi,d) T(theta,phi,d) TF(w) GaussIn(w)

f(t) = (power/sqrt(2 * pi * tau.^2)) * exp( -t^2 / (2 * tau.^2) ); % normallized gaussian input in time
fHan = matlabFunction(f(t));

F(w) = fourier(f,t,w); % the corresponding Fourier space function
FHan = matlabFunction(F(w));

R(theta) = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Pol(theta) = R(theta)*[1 0; 0 0]*R(theta);

BP(phi,d) = [exp(1i * ((dn/2)/c) * w * d)*exp(1i*phi) 0; 0 exp(1i * ((-dn/2)/c) * w * d)]; % Matrix of the actual crystal functions

T(theta,phi,d) = R(-theta)*BP(phi,d)*R(theta); % Rotated crystal to the right orientation


%% "Propagation"

GaussIn = [F(w);0]; % input pulse in frequency space

% Output of the stacker. There are currently four crystals. Also I don't
% call it T(w) right away so that the first element (along the polarizer
% direction) can be pulled out
TF  = Pol(thetaPols(4)) * T(thetaCrys(4),phiCrys(4),dCrys(4))...
    * Pol(thetaPols(4)) * T(thetaCrys(3),phiCrys(3),dCrys(3))...
    * Pol(thetaPols(4)) * T(thetaCrys(2),phiCrys(2),dCrys(2))...
    * Pol(thetaPols(4)) * T(thetaCrys(1),phiCrys(1),dCrys(1))...
    * GaussIn;
TF(w) = TF(1); % now we can make out function
TFHan = matlabFunction(TF(w));


Pulse(t) = ifourier(TF,w,t); % output pulse
if double(Pulse(0)) == 0
    PulseHan = @(t)0*t;
else
    PulseHan = matlabFunction(Pulse(t));
end




%% Plotting and analysis

dt = 0.01*ps;
times = -20*ps:dt:20*ps;

[FWHM{1},rise{1},~] = pulseParams(abs(fHan(times)).^2,times);
disp(['>>>> Input pulse has a FWHM of: ',num2str(FWHM{1}),...
    ' ps and a rise/fall time of: ',num2str(rise{1}),' ps'])

[FWHM{2},rise{2},~] = pulseParams(abs(PulseHan(times)).^2,times);
if length(rise{2}) >= 2
    disp(['>>>> Output pulse is a train with FWHM of: ',num2str(FWHM{2}(1)),...
        ' ps and a rise/fall time of: ',num2str(rise{2}(1)),' ps'])
else
    disp(['>>>> Output pulse has a FWHM of: ',num2str(FWHM{2}),...
        ' ps and a rise/fall time of: ',num2str(rise{2}),' ps'])
end

figure(figNum)
clf

subplot(2,1,1)
plot(times,abs(fHan(times)).^2,'LineWidth',3)
title('Input Pulse','FontSize',30);
xlabel('ps','FontSize',24);
text(min(times)*0.925,max(abs(fHan(times)).^2)*0.85,'Intital Parameters','FontWeight','bold');
text(min(times)*0.925,max(abs(fHan(times)).^2)*0.8,['Crystal Choice: ',crys]);
text(min(times)*0.925,max(abs(fHan(times)).^2)*0.75,['Wavelength: ',num2str(lambda),' um']);
text(min(times)*0.925,max(abs(fHan(times)).^2)*0.7,['Crystal Lengths: ',char(strjoin(string(dCrys),', ')),' mm']);
text(min(times)*0.925,max(abs(fHan(times)).^2)*0.65,['FWHM: ',num2str(FWHM{1}),' ps']);
text(min(times)*0.925,max(abs(fHan(times)).^2)*0.6,['Rise/Fall: ',num2str(rise{1}),' ps']);
fillFig(0,0)

subplot(2,1,2)
plot(times,abs(PulseHan(times)).^2,'LineWidth',3)
title('Output Pulse','FontSize',30);
xlabel('ps','FontSize',24);
text(min(times)*0.925,max(abs(PulseHan(times)).^2)*0.85,'Output Parameters','FontWeight','bold');
text(min(times)*0.925,max(abs(PulseHan(times)).^2)*0.8,['FWHM: ',num2str(FWHM{2}(1)),' ps']);
text(min(times)*0.925,max(abs(PulseHan(times)).^2)*0.75,['Rise/Fall: ',num2str(rise{2}(1)),' ps']);
fillFig(0,0)



%% Bonus functions

function [ne,no] = nonLinCrysChoice(crysType,lambda)

% Sellmeier equations

ne_BBO = sqrt( 2.31197 + 0.01184/(lambda^2 - 0.01607) - 0.00400 * lambda^2 );
no_BBO = sqrt( 2.67579 + 0.02099/(lambda^2 - 0.00470) - 0.00528 * lambda^2 );

ne_CAL = sqrt( 2.18438 + 0.0087309/(lambda^2 - 0.01018) - 0.0024411 * lambda^2 );
no_CAL = sqrt( 2.69705 + 0.0192064/(lambda^2 - 0.01820) - 0.0151624 * lambda^2 );

ne_YVO4 = sqrt(4.59905 + 0.110534/(lambda^2 - 0.04813) - 0.0122676 * lambda^2 );
no_YVO4 = sqrt(3.77834 + 0.069736/(lambda^2 - 0.04724) - 0.0108133 * lambda^2 );



switch crysType
    case 'BBO'
        ne = ne_BBO;
        no = no_BBO;
    case 'CAL'
        ne = ne_CAL;
        no = no_CAL;
    case 'YVO4'
        ne = ne_YVO4;
        no = no_YVO4;
end

end

function [FWHM,riseTime,fallTime] = pulseParams(pulseVec,times)

FWHM = pulsewidth(pulseVec,times);
riseTime = risetime(pulseVec,times);
fallTime = falltime(pulseVec,times);


end

function [] = fillFig(uFillW,uFillH)
    
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-uFillW;
ax_height = outerpos(4) - ti(2) - ti(4)-uFillH;
ax.Position = [left bottom ax_width ax_height];

end