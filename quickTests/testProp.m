%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                            Constants                              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear;
% clc;

% Units are in um throughout, including the initial beam properties
m = 10^3;
mm = 10^-3 * m;
cm = 10^-2 * m;
um = 10^-6 * m;
km = 10^3 * m;

lambda = 1.55 * um; % Wavelength of light

%% Beam Creation

beam = beamPropagation2D(lambda,4*cm,4*cm,2^12,'hex');
beamBak = beam.outputProperties2D('');

%% Prop and/or Plot

beam.inputProperties2D(beamBak);

z =  5.217 * m; % Distance to Propagate

tmp = beam.outputProperties2D('hex');

tmp.hex_BeamsOn = [1,1,1,1,1,1,1];
% tmp.hex_PhaseOffset = [0,1.4302,1.4568,1.2987,1.187,1.3087,1.4081]*pi;
tmp.hex_PhaseOffset = [pi,0,0,0,0,0,0];
tmp.hex_PhaseCurve = ones(1,7)*-1627;

% tmp.hex_PhaseCurve = repmat(2 * m,1,7);
% tmp.hex_DistBeams = 3.5 * mm;
% tmp.hex_AperX = tmp.hex_DistBeams - .05 * mm;
% tmp.hex_AperY = tmp.hex_AperX;
% beam.gen_nPlotPoints = 400;
beam.hex_InitialBeamDef2D(tmp);

showfigs = [1,1,1,1];

saveFlag = 0;
folderName = 'C:\Users\rlemons\Documents\Projects\ULM\PaperFigures\';
runName = 'combo1_V2_';
fileType = {'.svg','.fig'};


N = size(beam.field_fList,1)/2;
M = N-floor(N/3):N+floor(N/3)-1;

if showfigs(1)
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
end

if showfigs(2)
f(2) = figure(2);clf;
% ax2 = subplot(2,2,2);
% beamPhase = angle(beam.field_fList(M,M))+pi-pi/4;
% back = beamPhase == pi-pi/4;
% beamPhase(back) = 0;
% beamPhase(~back) = beamPhase(~back)+pi;
beamPhase = angle(beam.field_fList(M,M));

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
end


[beam.field_fList,beam.field_FList] = beam.forwardProp_FreeSpace2D(z);

N = size(beam.field_fList,1);
M = round(linspace(1,N,length(M)));

if showfigs(3)
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
end

if showfigs(4)
f(4) = figure(4);clf;
% ax4 = subplot(2,2,4);
% beamPhase = angle(beam.field_fList(M,M))+pi;
beamPhase = angle(beam.field_fList(M,M));

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
end

if saveFlag
    if showfigs(1)
    saveas(f(1),[folderName,runName,'nearInt',fileType{1}]);
    saveas(f(1),[folderName,runName,'nearInt',fileType{2}]);
    end
    if showfigs(2)
    saveas(f(2),[folderName,runName,'nearPhase',fileType{1}]);
    saveas(f(2),[folderName,runName,'nearPhase',fileType{2}]);
    end
    if showfigs(3)
    saveas(f(3),[folderName,runName,'farInt',fileType{1}]);
    saveas(f(3),[folderName,runName,'farInt',fileType{2}]);    
    end
    if showfigs(4)
    saveas(f(4),[folderName,runName,'farPhase',fileType{1}]);
    saveas(f(4),[folderName,runName,'farPhase',fileType{2}]);
    end
end