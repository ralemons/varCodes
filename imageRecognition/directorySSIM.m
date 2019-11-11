%% Load directory information
clear

% load in .mat files with directory information. (better git integration)
load('dirNames.mat');

%% Find same named images
% dirProps = dir(startDir);
% 
% % remove all directories from dir listing
% N = size(dirProps,1);
% for ii = N:-1:1
%     if dirProps(ii).isdir == 1
%         dirProps(ii) = [];
%     end
% end
% N = size(dirProps,1); % find new size of the file list
% 
% sameNames = zeros(1,N);
% 
% for ii = 1:N
%     if ~isempty(regexp(dirProps(ii).name,'\(\d*\)', 'once'))
%         sameNames(ii) = 1;
%     end
% end
% 
% indSame = find(sameNames)';
% indOther = zeros(length(indSame),1);
% 
% for ii = 1:length(indSame)
%     nameFind{ii,2} = dirProps(indSame(ii)).name;
%     tmp1 = extractBefore(dirProps(indSame(ii)).name,'(');
%     tmp2 = extractAfter(dirProps(indSame(ii)).name,')');
% %     if tmp1(end) == ' '
%         tmp1(end) = [];
% %     end
%     nameFind{ii,1} = [tmp1,tmp2];
% end
% 
% 
% tmp = struct2cell(dirProps)';
% tmp = tmp(:,1);
% 
% for ii = 1:length(nameFind)
%     if ~isempty(find(strcmp(tmp,nameFind{ii,1})))
%         indOther(ii) = find(strcmp(tmp,nameFind{ii,1}));
%     end
% end
% 
% 
% for ii = 1:length(indSame)
%    I1 = imread(fullfile(...
%        dirProps(indSame(ii)).folder,dirProps(indSame(ii)).name));
%    I1 = insertText(I1,[0 0; 0 0],dirProps(indSame(ii)).name,...
%        'FontSize',24);
%    I2 = imread(fullfile(...
%        dirProps(indOther(ii)).folder,dirProps(indOther(ii)).name));
%    I2 = insertText(I2,[0 0; 0 0],dirProps(indOther(ii)).name,...
%        'FontSize',24);
%    subplot(2,1,1);
%    imshow(I1)
%    subplot(2,1,2);
%    imshow(I2)
%    userIn = input('Actually Same?','s');
% end
% 
% for ii = 1:N
%     if sameNames(ii) == 1
%         movefile(fullfile(dirProps(ii).folder,dirProps(ii).name),moveDir)
%     end
% end

%% Load up the smaller list of files
clear sameNames dirPropsdirP

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
    [imgs{ii},cmap,~] = imread(fileName);
    if ~isempty(cmap)
        imgs{ii} = ind2rgb(imgs{ii},cmap);
        cmap = [];
    end
    if size(imgs{ii},3) == 3
        imgs{ii} = rgb2gray(imgs{ii}); % ssim needs two dimentional data
    end
    imgs{ii} = im2single(imresize(imgs{ii},[300 300]));
end

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
clf
ax = gca;

nComp = 8;

for ii = 1:N
    [vals,inds] = sort(ssimMat(ii,:));
    trash = vals == 1 | vals == 0;
    vals(trash) = [];
    inds(trash) = [];
    inds = flip(inds);
    
    inds = inds(1:nComp);
    
    
    imgsDisp{4} = 0;
    
    [imgsDisp{1},cmap,~] = imread(fullfile(dirProps(ii).folder,dirProps(ii).name));
    if ~isempty(cmap)
        imgsDisp{1} = ind2rgb(imgsDisp{1},cmap);
        cmap = [];
    end
    
    imgsDisp{1} = insertText(imgsDisp{1},...
        [0 0;0 round(size(imgsDisp{1},1).*0.08.*2)],dirProps(ii).name,...
        'AnchorPoint','LeftBottom',...
        'FontSize',min(200,round(size(imgsDisp{1},1).*0.08)));
    
    if dupeFiles(ii) == 1
        imgsDisp{1} = insertText(imgsDisp{1},...
        [round(size(imgsDisp{1},1)/2) round(size(imgsDisp{1},1)/2)],'DUPE',...
        'AnchorPoint','Center',...
        'FontSize',min(200,round(size(imgsDisp{1},1).*0.2)));
    end
    
    for jj = 1:nComp
        [imgsDisp{jj+1},cmap,~] =...
            imread(fullfile(dirProps(inds(jj)).folder,dirProps(inds(jj)).name));
        
        if ~isempty(cmap)
            imgsDisp{jj+1} = ind2rgb(imgsDisp{jj+1},cmap);
            cmap = [];
        end
        
        imgsDisp{jj+1} = insertText(imgsDisp{jj+1},...
            [0 0;0 round(size(imgsDisp{jj+1},1).*0.08.*2)],dirProps(inds(jj)).name,...
            'AnchorPoint','LeftBottom',...
            'FontSize',round(size(imgsDisp{jj+1},1).*0.08));
        
        if dupeFiles(inds(jj)) == 1
            imgsDisp{jj+1} = insertText(imgsDisp{jj+1},...
                [round(size(imgsDisp{jj+1},1)/2) round(size(imgsDisp{jj+1},1)/2)],'DUPE',...
                'AnchorPoint','Center',...
                'FontSize',min(200,round(size(imgsDisp{jj+1},1).*0.2)));
        end
        
    end
    
    montage(imgsDisp,'Parent',ax);
    
%     drawnow
%     pause(0.25);

    userAns = input('Too Similar?','s');
    if userAns == 'l'
        dupeFiles(ii) = 1;
    end
end

%% New names and/or directory
names = cell(N,1);

for ii = 1:N
    names{ii,1} = strrep(dirProps(ii).name,'png','gif');
end

%% Move Identical Files w/ Different Names
for ii = 1:N
    if dupeFiles(ii) == 1
        movefile(fullfile(otherDir,names{ii}),moveDir)
%         movefile(fullfile(dirProps(ii).folder,dirProps(ii).name),moveDir)
    end
end
disp('fin')
