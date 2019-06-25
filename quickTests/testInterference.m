%%% Two Beam Interference Test
%%% This is designed to test the applicability of using two propagations of
%%% the fields through the 2D beam prop as analogs for polarization. Then
%%% combining them to get the proper interference fringes.

%% Create Beams

clear

hBeam = beamPropagation2D(.00155,12,12,2^11,'hex');
vBeam = beamPropagation2D(.00155,12,12,2^11,'hex');

% hBeam = beamPropagation2D(.00155,12,12,2^11,'rect');
% vBeam = beamPropagation2D(.00155,12,12,2^11,'rect');

save('testVars.mat');
    

%% Load and set parameters

clear
load('testVars.mat');

% hBeam.field_Polar = [1 0];
% vBeam.field_Polar = [0 1];

hBeam.field_Polar = (1/sqrt(2)) * [1 -1i];
vBeam.field_Polar = (1/sqrt(2)) * [1 1i];


%% Propagate

[hBeam.field_fList,hBeam.field_FList] = hBeam.forwardProp_FreeSpace2D(5000);
[vBeam.field_fList,vBeam.field_FList] = vBeam.forwardProp_FreeSpace2D(5000);



%% Get Values

r1 = abs(vBeam.field_fList);
ph1 = angle(vBeam.field_fList);
po1 = vBeam.field_Polar;

r2 = abs(hBeam.field_fList);
ph2 = angle(hBeam.field_fList);
po2 = hBeam.field_Polar;

%% Plot

inter = r1.^2 + r2.^2 + 2*r1.*r2.*cos(ph1-ph2)*dot(po1,po2);

figure(1);
pcolor(inter);
shading interp;
title('Orthogonally Polarized');
shading interp;
figure(2);
pcolor(abs(hBeam.field_fList + vBeam.field_fList).^2);
shading interp;
title('Parallel Polarized')
shading interp;