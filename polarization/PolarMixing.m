%% Derive the Stokes parameters of a beam from the photos %%
%
% Faily basic code which takes seven images, slices them up into 
% macropixels (MPs), and then generates the stokes parameters for each MP.
%
% Assumes that the group of seven images is in a single folder alone and
% that they are .png's. Works best (at all) if the images are black and
% white to start with. If they are not then matlab converts them to BW with
% the built in weighting scheme. This scheme will not properly convert
% false color images to intesity properly.
%
% For best results you should hand crop the images so that any offsets sue
% to the optics moving between measurements is mitigated.
%

%% Load the folder
clear

% CHANGE THIS LINE
filePath = uigetdir('D:\ULM\Polarization\Pictures\'); % points to the general folder holding the images


%% Load the files

dataIN(:,:,1) = imread(fullfile(filePath,'Full-IR.png'));
dataIN(:,:,2) = imread(fullfile(filePath,'H-IR.png'));
dataIN(:,:,3) = imread(fullfile(filePath,'V-IR.png'));
dataIN(:,:,4) = imread(fullfile(filePath,'D-IR.png'));
dataIN(:,:,5) = imread(fullfile(filePath,'AD-IR.png'));
dataIN(:,:,6) = imread(fullfile(filePath,'L-IR.png'));
dataIN(:,:,7) = imread(fullfile(filePath,'R-IR.png'));

% dataIN(:,:,1) = rgb2gray(imread(fullfile(filePath,'Full-IR.png')));
% dataIN(:,:,2) = rgb2gray(imread(fullfile(filePath,'H-IR.png')));
% dataIN(:,:,3) = rgb2gray(imread(fullfile(filePath,'V-IR.png')));
% dataIN(:,:,4) = rgb2gray(imread(fullfile(filePath,'D-IR.png')));
% dataIN(:,:,5) = rgb2gray(imread(fullfile(filePath,'AD-IR.png')));
% dataIN(:,:,6) = rgb2gray(imread(fullfile(filePath,'L-IR.png')));
% dataIN(:,:,7) = rgb2gray(imread(fullfile(filePath,'R-IR.png')));

% dataIN(:,:,1) = imread(fullfile(filePath,'Full-IR.bmp'));
% dataIN(:,:,2) = imread(fullfile(filePath,'H-IR.bmp'));
% dataIN(:,:,3) = imread(fullfile(filePath,'V-IR.bmp'));
% dataIN(:,:,4) = imread(fullfile(filePath,'D-IR.bmp'));
% dataIN(:,:,5) = imread(fullfile(filePath,'AD-IR.bmp'));
% dataIN(:,:,6) = imread(fullfile(filePath,'L-IR.bmp'));
% dataIN(:,:,7) = imread(fullfile(filePath,'R-IR.bmp'));


%% Divide them into smaller groupings and Process

clearvars -except dataIN filePath

dataOUT = double(dataIN);

% for ii = 1:size(dataOUT,3)
%    
%     dataOUT(:,:,ii) = dataOUT(:,:,ii)./max(max(dataOUT(:,:,ii)));
%     
% end


numGrid.H = 30; % number of MPs in x
numGrid.V = 30; % number of MPs in y

% Find out number of pixels per MP in x and y
numPix.H = floor(size(dataOUT(:,:,1),2)/numGrid.H);
numPix.V = floor(size(dataOUT(:,:,1),1)/numGrid.V);

% Based on number of MP and number of pixels per MP, calculates the window
% size. This will always be less than or equal to the original image size
sizeWind.H = numPix.H*numGrid.H;
sizeWind.V = numPix.V*numGrid.V;

% Start position of this smaller window
posWind.H = size(dataOUT(:,:,1),2)/2-sizeWind.H/2;
posWind.V = size(dataOUT(:,:,1),1)/2-sizeWind.V/2;

% Calculates which rows and columns to elliminate
del.H = [1:posWind.H 1+sizeWind.H+posWind.H:size(dataOUT(:,:,1),2)];
del.V = [1:posWind.V 1+sizeWind.V+posWind.V:size(dataOUT(:,:,1),1)];

% Handles odd cases where the windows is equal to the image size in x and y
if ~isempty(del.H) && ~isempty(del.V)
    dataOUT(del.V,:,:) = [];
    dataOUT(:,del.H,:) = [];
elseif isempty(del.H) && ~isempty(del.V)
    dataOUT(del.V,:,:) = [];
elseif ~isempty(del.H) && isempty(del.V)
    dataOUT(:,del.H,:) = [];
end

% Generate the indexing numbers in x and y for each MP
for ii = 1:numGrid.H
    grids.H(:,ii) = 1+(numPix.H*(ii-1)):numPix.H+(numPix.H*(ii-1));
end
for ii = 1:numGrid.V
    grids.V(:,ii) = 1+(numPix.V*(ii-1)):numPix.V+(numPix.V*(ii-1));
end


% Makes intesity of each MP equal to mean of all pixels inside it
for kk = 1:size(dataOUT,3)
    for ii = 1:numGrid.V
        for jj = 1:numGrid.H
            
            dataOUT(grids.V(:,ii),grids.H(:,jj),kk) =...
                mean2(dataOUT(grids.V(:,ii),grids.H(:,jj),kk));
            
        end
    end
end


% This is the stokes parameters but of the image size ie. 1024x1280
S = zeros(sizeWind.V,sizeWind.H,4);

% This is the stokes parameters but of the number of MP in size ie. 30x30
S_num = zeros(numGrid.V,numGrid.H,4);

% Element wise subtraction of MP
for ii = 1:numGrid.V
    for jj = 1:numGrid.H
        
        p0 = dataOUT(grids.V(:,ii),grids.H(:,jj),1);
        p1 = (dataOUT(grids.V(:,ii),grids.H(:,jj),2) - dataOUT(grids.V(:,ii),grids.H(:,jj),3));
        p2 = (dataOUT(grids.V(:,ii),grids.H(:,jj),4) - dataOUT(grids.V(:,ii),grids.H(:,jj),5));
        p3 = (dataOUT(grids.V(:,ii),grids.H(:,jj),6) - dataOUT(grids.V(:,ii),grids.H(:,jj),7));
        
        
        S(grids.V(:,ii),grids.H(:,jj),1) = p0; % should be p0/p0 but this gives the orignal image and p0 isn't used
        S(grids.V(:,ii),grids.H(:,jj),2) = p1./p0;
        S(grids.V(:,ii),grids.H(:,jj),3) = p2./p0;
        S(grids.V(:,ii),grids.H(:,jj),4) = p3./p0;
        
        S_num(ii,jj,1) = unique(p0); % same as with S
        S_num(ii,jj,2) = unique(p1./p0);
        S_num(ii,jj,3) = unique(p2./p0);
        S_num(ii,jj,4) = unique(p3./p0);
                
        
    end
end

% Create the sum of the squares to get visual of error. Farther from 1 in 
% either direction is larger error
S(:,:,5) = sum( S(:,:,2:4).^2 , 3);
S_num(:,:,5) = sum( S_num(:,:,2:4).^2 , 3);

% tmp = S(:,:,5);
% tmp_num = S_num(:,:,5);
% 
% range = .2;
% tmp(tmp<(1+range) & tmp>(1-range)) = 100;
% tmp_num(tmp_num<(1+range) & tmp_num>(1-range)) = 100;
% 
% S(:,:,5) = tmp;
% S_num(:,:,5) = tmp_num;

% Plot the stokes parameter
figure(1)
for ii = 1:3
    subplot(2,2,ii)
    imagesc(S_num(:,:,ii+1))
%     caxis([-1,1]) % uncomment to enforce color limits. Another error
    colorbar
    axis square
    fillFig(0,0)
    title(sprintf('S_%i',ii),'Interpreter','tex')
end

% Plot the sum of the squares
subplot(2,2,4)
imagesc(S_num(:,:,5))
colorbar
axis square
fillFig(0,0)
title('\Sigma S_i^2','Interpreter','tex')



%% Write out the data to csv
% 
% nRows = numel(S_num(:,:,1));
% 
% % The out var is 6 diminsioned to work with the poincare.m function from
% % the matlab file exchange.
% out = zeros(nRows,6);
% out(:,1) = 1:nRows;
% for ii = 1:4
%     tmp = S_num(:,:,ii);
%     out(:,ii+1) = tmp(:);   
% end
% 
% % Auto save to a file name based on in the input folder
% tmp = split(filePath,'\');
% tmp = tmp{end};
% fileOutPath = fullfile('D:\ULM\Polarization\Results\',[tmp '.csv']);
% 
% writematrix(out,fileOutPath); % You can then load this directly by the .m
% 

%% Plot the polarization ellipses over the full image of the beam

warning('off','MATLAB:plot:IgnoreImaginaryXYPart');

clear S_polar

% Take only S1, S2, and S3
S_polar = S_num(:,:,2:4);

fullIm = dataIN(:,:,1);

div = 2; % Divides number of MP by this to choose how many MP to plot
scale = 75; % How big each arrrow is
arrowMove = 1200; % How far the arrow is along the circle
fignum = 2; % What the number of the figure you create is
chooseROI = 1; % Pulls up selectable ROI if 1, plot full image if 0 

polarEllipsePlot(S_polar,fullIm,div,scale,arrowMove,fignum,chooseROI);

