%% This script is designed to be used to simualate steering the beam
%
% Right now this script is functional in the sense that is simulates
% exactly what we want it to. In the future, expanding this code to either
% output phase offsets for the ULM or directly driving the ULM with the
% offset is the goal.
%
%% Create the beam object and other constants

clear

% Unit Definitions & Propagation Constants
m = 10^3; % If m is not 1, then it is in the reciprocal unit space. ie. m = 10^6 means units of um (10^-6);
cm = m * 10^-2;
mm = m * 10^-3;
um = m * 10^-6;
nm = m * 10^-9;

% Distance to propagate. Everything was built and hard coded around 2.5m
zProp = 25 * m;
lambda = 1.55 * um; % Wavelength of laser. ULM is 1.55

xSize = 6 * cm;
ySize = 6 * cm;


%% Beam Definition


% Beam object to be used
beam = beamPropagation2D(lambda,xSize,ySize,2^7,'hex');
params = beam.outputProperties2D('hex');

% Save the Object
fileLoc = 'D:\ULM\Beam Steering';
fileName = 'S_in_SLAC_Beam';
save(fullfile(fileLoc,fileName),'beam');


%% Load the Object (So you don't have to create a beam each time)

fileLoc = 'D:\ULM\Beam Steering';
fileName = 'S_in_SLAC_Beam';
load(fullfile(fileLoc,fileName));


%% Create offsets for steering

% Location of folder containing alphabet files
alphaLoc = 'D:\ULM\Beam Steering\Alphabet\';

word = 'W'; % Word to generate
nP = 50; % Number of phase offsets per letter
phaseOffs = zeros(nP,7,length(word)); % Phase offsets for each beam
x = zeros(nP,length(word));
y = zeros(nP,length(word));

% These are the phase offsets that generate max shifts in H or V
% positive direction. By making each offset a linear combination of a
% weighting of hPCol and vPCol we get the offset needed for any shift
hPCol = [0, 0, 2*pi/3, 2*pi/3, 0, -2*pi/3, -2*pi/3];
% hPCol = [0, 0, 0, 0, 0, 0, 0];
vPCol = [0, -2*pi/3, -pi/3, pi/3, 2*pi/3, pi/3, -pi/3];
% vPCol = [0, 0, 0, 0, 0, 0, 0];

for ii = 1:length(word)
    
    % Get the letter file
    alphaFile = fullfile([alphaLoc,word(ii),'.png']);
    [~,~,I] = imread(alphaFile);
    I = imbinarize(I);
    
    % Generate path of the outside boundary
    tmp = bwboundaries(flip(I',2));
    mask = tmp{1};
    
    % Takes only nP points in the path
    mask = mask(1:length(mask)/nP:end,:);
    
    % Generate x and y lists of points, center at 0 scaled from -1 to 1
    x(:,ii) = mask(:,1)/(size(I,1)/2);
    x(:,ii) = x(:,ii) - 1;
    y(:,ii) = mask(:,2)/(size(I,2)/2);
    y(:,ii) = y(:,ii) - 1;
    
    % Actually builds the offset
    phaseOffs(:,:,ii) = x(:,ii) .* hPCol + y(:,ii) .* vPCol;
    
end


%% Propagation, plotting, and giffing


% Very hardcoded zoomed window for display. Makes it look nicer
sx = .13; % zero to 1
sy = .14; % zero to 1
N = beam.grid_npts;
plotSizeX = (N/2 + 1)-floor(sx*N/2):(N/2 + 1)+floor(sx*N/2);
plotSizeY = (N/2 + 1)-floor(sy*N/2):(N/2 + 1)+floor(sy*N/2);


for kk = 1:1 % Loop to repeat the movement 
    
    for ii = 1:size(phaseOffs,1) % Loop over the phase offsets (nP)
        
        for jj = 1:size(phaseOffs,3) % Loop over letters
            
            % Set the beam as having given phase shift
            params.hex_PhaseOffset = phaseOffs(ii,:,jj);
            beam.hex_InitialBeamDef2D(params);
            
            % Propagate the beam
            a = beam.forwardProp_FreeSpace2D(zProp);
            

            
            % Plot the prop'd beam and an overlay of the letter it is
            % supposed to be
            subplot(1,size(phaseOffs,3),jj)
            imagesc([-1 1],[-1 1],flipud(abs(a(plotSizeY,plotSizeX)).^2))
            axis square
            hold on
            scatter( x(:,jj),-y(:,jj),...
                50,[140,21,21]/255,'filled')
            hold off
            set(gca,'visible','off')
            
%             pause(0.1) %pause to see it evolve

%             imagesc(abs(a).^2) %different plotting for whole space
%             axis square


        end
        
        drawnow % Makes it into a movie instead of only plotting last one
        
    end
    
end


