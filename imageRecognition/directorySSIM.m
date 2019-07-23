%% Load directory and image information
clear

% load in .mat files with directory information. (better git integration)
load('dirNames.mat');
dirProps = dir(startDir);

% remove all directories from dir listing
N = size(dirProps,1);
for ii = N:-1:1
    if dirProps(ii).isdir == 1
        dirProps(ii) = [];
    end
end
N = size(dirProps,1); % find new size of the file list

[imgs{1:N}] = deal(zeros(300)); % pre-allocate for thumbnail cell array
for ii = 1:N
    fileName = fullfile(dirProps(ii).folder,dirProps(ii).name);
    imgs{ii} = imresize(imread(fileName),[300 300]);
    if size(imgs{ii},3) == 3
        imgs{ii} = rgb2gray(imgs{ii}); % ssim needs two dimentional data
    end
end

% imgs = repmat(imgs,1,10);
% N = size(imgs,2);

%% Compare the images
ssimMat = zeros(N);
dupeFiles = zeros(N,1);


for jj = 1:N % Iterate through each image
    for ii = jj+1:N % Compare only to images after it
        ssimMat(jj,ii) = slimSSIM(imgs{jj},imgs{ii});
        if ssimMat(jj,ii) == 1
            dupeFiles(ii) = 1; % record if there is a truly identical image
        end
    end
end
ssimMat = ssimMat + ssimMat';

%% Display very close files for comparison

figure(1)
ax = gca;

for ii = 1:N
    [vals,inds] = sort(ssimMat(ii,:));
    trash = vals == 1 | vals == 0;
    vals(trash) = [];
    inds(trash) = [];
    inds = flip(inds);
    
    imgsDisp{4} = 0;
    
    imgsDisp{1} = imread(fullfile(dirProps(ii).folder,dirProps(ii).name));
    imgsDisp{1} = insertText(imgsDisp{1},...
        [0 0;0 round(size(imgsDisp{1},1).*0.032.*2)],dirProps(ii).name,...
        'AnchorPoint','LeftBottom',...
        'FontSize',round(size(imgsDisp{1},1).*0.032));
    for jj = 1:3
        imgsDisp{jj+1} =...
            imread(fullfile(dirProps(inds(jj)).folder,dirProps(inds(jj)).name));
        imgsDisp{jj+1} = insertText(imgsDisp{jj+1},...
            [0 0;0 round(size(imgsDisp{jj+1},1).*0.032.*2)],dirProps(inds(jj)).name,...
            'AnchorPoint','LeftBottom',...
            'FontSize',round(size(imgsDisp{jj+1},1).*0.032));
    end
    
    montage(imgsDisp,'Parent',ax);
    userAns = input('Too Similar?','s');
    
    if userAns == 'l'
        disp(1)
    end
end

%% Move Identical Files
for ii = 1:N
    if dupeFiles(ii) == 1
        movefile(fullfile(dirProps(ii).folder,dirProps(ii).name),moveDir)
    end
end

