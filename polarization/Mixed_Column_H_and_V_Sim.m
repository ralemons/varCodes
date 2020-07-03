%%% Two Beam Interference Test
%%% This is designed to test the applicability of using two propagations of
%%% the fields through the 2D beam prop as analogs for polarization. Then
%%% combining them to get the proper interference fringes.

%% Create Beams

clear

% Units are in um throughout, including the initial beam properties
m = 10^3;
mm = 10^-3 * m;
cm = 10^-2 * m;
um = 10^-6 * m;
km = 10^3 * m;

lambda = 1.5 * um; % Wavelength of light

beam = beamPropagation2D(lambda,8*cm,8*cm,2^11,'hex');
hBeam = beamPropagation2D();
vBeam = beamPropagation2D();
beamBak = beam.outputProperties2D('');


%% Load and set parameters

hBeam.inputProperties2D(beamBak);
vBeam.inputProperties2D(beamBak);

hBeam.field_Polar = [1 0];
% vBeam.field_Polar = [0 1];
vBeam.field_Polar = [1 0];

% hBeam.field_Polar = (1/sqrt(2)) * [1 -1i];
% vBeam.field_Polar = (1/sqrt(2)) * [1 1i];

z =  20* m; % Distance to Propagate
showFigs = [0,0,1,1];

beam.inputProperties2D(beamBak);
tmpH = beam.outputProperties2D('hex');
tmpV = beam.outputProperties2D('hex');

tmpH.hex_BeamsOn = [0,0,1,1,0,1,1];
tmpH.hex_PhaseOffset = [0,0,0,0,0,0,0] + pi/4;
% tmpH.hex_PhaseOffset = [0,0,-pi/2,-pi/2,0,pi/2,pi/2] + pi/4;

tmpH.hex_PhaseCurve = repmat(2 * m,1,7);
tmpH.hex_DistBeams = 3.5 * mm;
tmpH.hex_AperX = tmpH.hex_DistBeams - .05 * mm;
tmpH.hex_AperY = tmpH.hex_AperX;
hBeam.gen_nPlotPoints = 400;
hBeam.hex_InitialBeamDef2D(tmpH);


tmpV.hex_BeamsOn = [1,1,0,0,1,0,0];
% tmpV.hex_PhaseOffset = [0,0,0,0,0,0,0] + pi/4;
tmpV.hex_PhaseOffset = [-pi,0,0,0,0,0,0] + pi/4;

tmpV.hex_PhaseCurve = repmat(2 * m,1,7);
tmpV.hex_DistBeams = 3.5 * mm;
tmpV.hex_AperX = tmpV.hex_DistBeams - .05 * mm;
tmpV.hex_AperY = tmpV.hex_AperX;
vBeam.gen_nPlotPoints = 400;
vBeam.hex_InitialBeamDef2D(tmpV);


%%%% Propagate %%%%

[hBeam.field_fList,hBeam.field_FList] = hBeam.forwardProp_FreeSpace2D(z);
[vBeam.field_fList,vBeam.field_FList] = vBeam.forwardProp_FreeSpace2D(z);

%%%% Get Values %%%%

r1 = abs(hBeam.field_fList);
ph1 = angle(hBeam.field_fList);
po1 = hBeam.field_Polar;

r2 = abs(vBeam.field_fList);
ph2 = angle(vBeam.field_fList);
po2 = vBeam.field_Polar;

%%%% Plot %%%%

inter = r1.^2 + r2.^2;

N = size(beam.field_fList,1)/2;
M = N-floor(N/3):N+floor(N/3)-1;

if showFigs(1) == 1
    f(1) = figure(1);
    drawnow
    figure(1);
    clf
    imagesc(abs(hBeam.field_fList + vBeam.field_fList).^2);
    set(gca,'YDir','normal')
    axis square
    title('Initial Intesity')
    fillFig(0,0)
end


beamPhase = angle(hBeam.field_fList + vBeam.field_fList)+pi-pi/4;
back = beamPhase == pi-pi/4;
beamPhase(back) = 0;
beamPhase(~back) = beamPhase(~back)+pi;

if showFigs(2) == 1
    f(2) = figure(2);
    drawnow
    figure(2);
    clf
    imagesc(beamPhase);
    set(gca,'YDir','normal')
    axis square
    title('Initial Phase')
    colorTicks = linspace(min(min(beamPhase)),max(max(beamPhase)),7);
    if sum(colorTicks) == 0
        colorTicks = linspace(0,2*pi,7);
    end
    colorbar('eastoutside',...
        'Ticks',colorTicks,...
        'TickLabels',{'$0$','$\pi/3$','$2\pi/3$','$\pi$','$4\pi/3$','$5\pi/3$','$2\pi$'},...
        'TickLabelInterpreter','latex')
    fillFig(0.005,0)
end

N = size(beam.field_fList,1);
M = round(linspace(1,N,length(M)));


if showFigs(3) == 1
    f(3) = figure(3);
    drawnow
    figure(3);
    clf
    imagesc(inter)
    set(gca,'YDir','normal')
    axis square
    title('Far Field, Orthogonally Polarized, Intensity')
    fillFig(0.015,0.0)
end

if showFigs(4) == 1
    f(4) = figure(4);
    drawnow
    figure(4);
    clf
    imagesc(abs(hBeam.field_fList + vBeam.field_fList).^2);
    set(gca,'YDir','normal')
    axis square
    title('Far Field, Parallelly Polarized, Intensity')
    fillFig(0,0)
end
