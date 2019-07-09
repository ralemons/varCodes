%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                            Constants                              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

% Units are in um throughout, including the initial beam properties
m = 10^3;
mm = 10^-3 * m;
cm = 10^-2 * m;
um = 10^-6 * m;
km = 10^3 * m;

z =  25 * m; % Distance to Propagate
lambda = 1.5 * um; % Wavelength of light

%% Beam Creation

beam = beamPropagation2D(lambda,3*cm,3*cm,2^10,'hex');
beamBak = beam.outputProperties2D('');

%% Prop and/or Plot

beam.inputProperties2D(beamBak);

tmp = beam.outputProperties2D('hex');
tmp.hex_BeamsOn = [1,1,1,1,1,1,1];
% tmp.hex_AmpBeams = [1,.7,.7,1,1,1,.7];
% tmp.hex_PhaseOffset = [0,0,0,0,pi,pi,pi]+pi/6;
% tmp.hex_PhaseOffset = [0,pi,0,pi,0,pi,0]+pi/6;
% tmp.hex_PhaseOffset = [0,0,0,0,pi,0,0]+pi/6;
% tmp.hex_PhaseCurve = repmat(-1.5*m,1,7);
beam.hex_InitialBeamDef2D(tmp);

ax1 = subplot(2,2,1);
beam.plotField2D(beam.field_fList,'abs')
ax2 = subplot(2,2,2);
beam.plotField2D(beam.field_fList,'angle')


[beam.field_fList,beam.field_FList] = beam.forwardProp_FreeSpace2D(z);
% [beam.field_fList,beam.field_FList] = beam.backwardProp_FreeSpace2D(2*z/6);


ax3 = subplot(2,2,3);
beam.plotField2D(beam.field_fList,'abs')
ax4 = subplot(2,2,4);
beam.plotField2D(beam.field_fList,'angle')


% folder = 'D:\ULM\Beam Shaping\Sim Data\';
% imSave1 = getframe(ax1);
% imSave2 = getframe(ax2);
% imSave3 = getframe(ax3);
% imSave4 = getframe(ax4);
% imSave1 = imSave1.cdata;
% imSave2 = imSave2.cdata;
% imSave3 = imSave3.cdata;
% imSave4 = imSave4.cdata;
% imwrite(imSave1,[folder,'zeroPI_preProp.png']);
% imwrite(imSave2,[folder,'zeroPIPhase_preProp.png']);
% imwrite(imSave3,[folder,'zeroPI_postProp.png']);
% imwrite(imSave4,[folder,'zeroPIPhase_postProp.png']);
