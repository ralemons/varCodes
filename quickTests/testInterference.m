%%% Two Beam Interference Test
%%% This is designed to test the applicability of using two propagations of
%%% the fields through the 2D beam prop as analogs for polarization. Then
%%% combining them to get the proper interference fringes.

%% Create Beams

% Units are in um throughout, including the initial beam properties
m = 10^3;
mm = 10^-3 * m;
cm = 10^-2 * m;
um = 10^-6 * m;
km = 10^3 * m;

lambda = 1.5 * um; % Wavelength of light

if ~exist('hBeam','var')
        
    hBeam = beamPropagation2D(lambda,6*cm,6*cm,2^10,'hex');
    vBeam = beamPropagation2D(lambda,6*cm,6*cm,2^10,'hex');
    
    % hBeam = beamPropagation2D(.00155,12,12,2^10,'rect');
    % vBeam = beamPropagation2D(.00155,12,12,2^10,'rect');
    
    save('C:\Users\rlemons\Documents\GitHub\varCodes\savedVars\testVars.mat');
end


%% Load and set parameters

clear
load('C:\Users\rlemons\Documents\GitHub\varCodes\savedVars\testVars.mat');

hBeam.field_Polar = [1 0];
vBeam.field_Polar = [0 1];

% hBeam.field_Polar = (1/sqrt(2)) * [1 -1i];
% vBeam.field_Polar = (1/sqrt(2)) * [1 1i];


%% Propagate

[f1,F1] = hBeam.forwardProp_FreeSpace2D(10*m);
[f2,F2] = vBeam.forwardProp_FreeSpace2D(10*m);



%% Get Values

r1 = abs(f1);
ph1 = angle(f1);
po1 = hBeam.field_Polar;

r2 = abs(f2);
ph2 = angle(f2);
po2 = vBeam.field_Polar;

%% Plot

inter = r1.^2 + r2.^2 + 2*r1.*r2.*cos(ph1-ph2)*dot(po1,po2);

f(1) = figure(1);
pcolor(abs(hBeam.field_fList + vBeam.field_fList).^2);
shading interp;
axis square
title('Initial Intesity')

f(2) = figure(2);
pcolor(angle(hBeam.field_fList + vBeam.field_fList));
shading interp;
axis square
title('Initial Phase')

f(3) = figure(3);
pcolor(inter);
shading interp;
axis square
title('Far Field, Orthogonally Polarized, Intensity');

f(4) = figure(4);
pcolor(abs(f1 + f2).^2);
shading interp;
axis square
title('Far Field, Parallelly Polarized, Intensity')

