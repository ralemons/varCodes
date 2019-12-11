%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                            Constants                              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% Units are in um throughout, including the initial beam properties
m = 10^3;
mm = 10^-3 * m;
cm = 10^-2 * m;
um = 10^-6 * m;
km = 10^3 * m;

z =  10 * m; % Distance to Propagate
lambda = 1.5 * um; % Wavelength of light

%% Beam Creation

beam = beamPropagation2D(lambda,5*cm,5*cm,2^12,'hex');
beamBak = beam.outputProperties2D('');

%% Prop and/or Plot

beam.inputProperties2D(beamBak);

tmp = beam.outputProperties2D('hex');
tmp.hex_BeamsOn = [1,1,1,1,1,1,1];
tmp.hex_PhaseOffset = [0,0,0,0,0,0,pi]+pi/4;
beam.gen_nPlotPoints = 400;
beam.hex_InitialBeamDef2D(tmp);

folderName = '~/Documents/research/SLAC/Fall19/ULM/paperFigures/';
runName = 'combo8_';
fileType = {'.svg','.fig'};
saveFlag = 1;

N = size(beam.field_fList,1)/2;
M = N-floor(N/3):N+floor(N/3)-1;

f(1) = figure(1);clf;
% ax1 = subplot(2,2,1);
beamInt = abs(beam.field_fList(M,M)).^2;
beamInt = beamInt/max(max(beamInt));

imagesc(beamInt)
axis off
daspect([1 1 1])
colorbar('eastoutside',...
    'Ticks',linspace(0,1,7),...
    'TickLabels',{'$0$','$0.17$','$0.33$','$0.50$','$0.67$','$0.83$','$1.0$'},...
    'TickLabelInterpreter','latex')
caxis([0 1])
% ax1.FontSize = 24;
set(gca,'FontSize',40);

fillFig(0.1,0)


f(2) = figure(2);clf;
% ax2 = subplot(2,2,2);
beamPhase = angle(beam.field_fList(M,M))+pi-pi/4;
back = beamPhase == pi-pi/4;
beamPhase(back) = 0;
beamPhase(~back) = beamPhase(~back)+pi;

imagesc(beamPhase)
axis off
daspect([1 1 1])
colorbar('eastoutside',...
    'Ticks',linspace(min(min(beamPhase)),max(max(beamPhase)),7),...
    'TickLabels',{'$0$','$\pi/3$','$2\pi/3$','$\pi$','$4\pi/3$','$5\pi/3$','$2\pi$'},...
    'TickLabelInterpreter','latex')
% ax2.FontSize = 24;
set(gca,'FontSize',40);

fillFig(0.1,0)


[beam.field_fList,beam.field_FList] = beam.forwardProp_FreeSpace2D(z);

N = size(beam.field_fList,1);
M = round(linspace(1,N,length(M)));


f(3) = figure(3);clf;
% ax3 = subplot(2,2,3);
beamInt = abs(beam.field_fList(M,M)).^2;
beamInt = beamInt/max(max(beamInt));

imagesc(beamInt)
axis off
daspect([1 1 1])
colorbar('eastoutside',...
    'Ticks',linspace(0,1,7),...
    'TickLabels',{'$0$','$0.17$','$0.33$','$0.50$','$0.67$','$0.83$','$1.0$'},...
    'TickLabelInterpreter','latex')
caxis([0 1])
% ax3.FontSize = 24;
set(gca,'FontSize',40);

fillFig(0.1,0)

f(4) = figure(4);clf;
% ax4 = subplot(2,2,4);
beamPhase = angle(beam.field_fList(M,M))+pi;

imagesc(beamPhase)
axis off
daspect([1 1 1])
colorbar('eastoutside',...
    'Ticks',linspace(min(min(beamPhase)),max(max(beamPhase)),7),...
    'TickLabels',{'$0$','$\pi/3$','$2\pi/3$','$\pi$','$4\pi/3$','$5\pi/3$','$2\pi$'},...
    'TickLabelInterpreter','latex')
% ax4.FontSize = 24;
set(gca,'FontSize',40);

fillFig(0.1,0)

if saveFlag
    saveas(f(1),[folderName,runName,'nearInt',fileType{1}]);
    saveas(f(2),[folderName,runName,'nearPhase',fileType{1}]);
    saveas(f(3),[folderName,runName,'farInt',fileType{1}]);
    saveas(f(4),[folderName,runName,'farPhase',fileType{1}]);
    saveas(f(1),[folderName,runName,'nearInt',fileType{2}]);
    saveas(f(2),[folderName,runName,'nearPhase',fileType{2}]);
    saveas(f(3),[folderName,runName,'farInt',fileType{2}]);
    saveas(f(4),[folderName,runName,'farPhase',fileType{2}]);
end