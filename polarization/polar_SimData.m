clear

% CHANGE THESE TWO LINES
rgbImage = imread('C:Projects\ULM\tmp.jpg'); % This line should point to any image
saveLoc = 'D:\ULM\Polarization\Pictures\Simulated_04262019\'; % This line should point to where you want to save the files

% Create full image
bwImage = rgb2gray(rgbImage);
imwrite(bwImage,[saveLoc,'Full-IR.png']);

% Hand built polarization weightings for each portion of the picture. They
% are built from back-solving the stokes parameter equations for a desired
% ellipse.
coeffs(:,1) = [1,0,0,0,0,0,0,0,0.065];      %H
coeffs(:,2) = [0,1,0,0,0,0,0,0,0.935];      %V
coeffs(:,3) = [0,0,1,0,0,0,0.935,0.065,0];  %D
coeffs(:,4) = [0,0,0,1,0,0,0.065,0.935,0];  %AD
coeffs(:,5) = [0,0,0,0,1,0,0.75,0.25,0.25]; %L
coeffs(:,6) = [0,0,0,0,0,1,0.25,0.75,0.75]; %R

% File names to work with PolarMixing.m
names = {'H-IR.png','V-IR.png','D-IR.png','AD-IR.png','L-IR.png','R-IR.png'};

% Build the remaining images
[Y,X] = size(bwImage);
for ii = 1:6
    
    
    tmp = bwImage;
    counter = 1; % Counts which region you are in
    
    for yy = 1:3
        for xx = 1:3
            
            % Generate weighting of image in desired region
            tmp( (((yy-1)/3)*Y)+1:(yy/3)*Y,(((xx-1)/3)*X)+1:(xx/3)*X ) =...
                tmp( (((yy-1)/3)*Y)+1:(yy/3)*Y,(((xx-1)/3)*X)+1:(xx/3)*X ) .* coeffs(counter,ii);
            
            counter = counter + 1;
            
        end
    end
    
    % Write out each image based on the 
    imwrite( tmp,fullfile(saveLoc,names{ii}) );
    
end
