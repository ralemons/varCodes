%% Parse the data generated and saved by the GA

clear
%% Load the data
if ~exist('data','var')
    data = load('C:\Users\rlemons\Documents\GitHub\beamProp\finalBeams.mat');
    [beam_orig,theHerd] = reconBeams(data);
end


genTime = mean([data.sol(:).timeG]);
runTime = mean([data.sol(:).timeZ]);
totTime = sum([data.sol(:).timeZ]);

fprintf(' \n')
disp(['>> Problem solved with ', upper(data.beamType), ' beam type',...
    ' for ', num2str(data.numGens), ' generations',...
    ' with ', num2str(data.numHerd), ' individuals.'])
disp(['>> The variable properties were: ', strjoin(data.propTest(:,1),', '),'.'])
disp(['>> Each generation took ~', num2str(round(genTime)),' s.'])
disp(['>> Each run took ~', num2str(round(runTime/60)),' min.'])
disp(['>> The total time taken was ', num2str(totTime/60/60),' h.'])
fprintf(' \n')

%% Concatinate error and sort by best runs

fitList = reshape(...
    cell2mat(...
    arrayfun(@(S)S(:).fitness(:,:),data.sol,'UniformOutput',false)...
    )...
    ,[],2,data.numRuns);
% fitList = fitList/max(max(max(fitList)));
fitList = fliplr(permute(permute(fitList,[1 3 2]),[2,1,3]));

[~,bestRuns(:,1)] = sortrows(fitList(:,:,1),'descend');
[~,bestRuns(:,2)] = sortrows(fitList(:,:,2),'descend');

%% Plot the solution data

f(1) = figure(1);
set(f,'Name','Beam to solve for');


if isa(beam_orig,'beamPropagation2D')
    subplot(1,2,1)
    beam_orig.plotField2D('abs');
    
    subplot(1,2,2)
    beam_orig.plotField2D('angle');
    
    minVal = min(min(abs(beam_orig.field_fList).^2));
    maxVal = max(max(abs(beam_orig.field_fList).^2));
% elseif strcmpi(data.solveType,'matrix')
    
else
    imagesc(beam_orig);
    axis square off
    
    minVal = min(min(beam_orig));
    maxVal = max(max(beam_orig));
end

%% Plot the Error lists

figure(50)
clf
clear lineColors


% Set Params
colorScale = 1.3;
colorBack = {'k','w'};
fontSize = [30,15];
smallVals = 10;



lineColors(:,:,1) = colormap(parula(data.numRuns*colorScale));
lineColors(:,:,2) = colormap(gray(data.numRuns*colorScale));
lineColors(1:end-data.numRuns,:,:) = []; lineColors = flipud(lineColors);

% lineColors(:,:,1) = flipud(colormap(parula(data.numRuns*colorScale)));
% lineColors(1:end-data.numRuns,:,:) = [];


hold on
for jj = 1:2
    for ii = data.numRuns:-1:1
        if ii == 1 || ii == data.numRuns
            plot(data.numGens:-1:1,fitList(bestRuns(ii,1),:,jj),...
                'color',lineColors(ii,:,jj),...
                'LineWidth',2)
        else
            plot(data.numGens:-1:1,fitList(bestRuns(ii,1),:,jj),...
                'color',lineColors(ii,:,jj))
        end
    end
end
hold off

title('Fitness Evolution of Each Generation');
set(gca,'FontSize',fontSize(1),'Color',colorBack{1});

xlim([1 data.numGens])


fillFig(0,0)

axes('Position' ,[.2 .6 .3 .3]);
box on

set(gca,'FontSize',fontSize(2),'Color',colorBack{1},...
    'XColor',colorBack{2},'YColor',colorBack{2});


xlim([1 smallVals])

hold on
for jj = 1:2
    for ii = data.numRuns:-1:1
        if ii == 1 || ii == data.numRuns
            plot(smallVals:-1:1,...
                fitList(bestRuns(ii,1),end-smallVals+1:end,jj),...
                'color',lineColors(ii,:,jj),...
                'LineWidth',2)
        else
            plot(smallVals:-1:1,...
                fitList(bestRuns(ii,1),end-smallVals+1:end,jj),...
                'color',lineColors(ii,:,jj))
        end
    end
end
hold off

%% All Sols Plot

clear plots I 
figure(99)

p = numSubPlots(data.numRuns);

for ii = 1:data.numRuns
    subplot(p(1),p(2),ii)
    
    plots(:,:,1) = abs(theHerd(ii).field_fList).^2;
    plots(:,:,2) = angle(theHerd(ii).field_fList) + pi;
    plots(:,:,2) = max(max(plots(:,:,1))) * (plots(:,:,2)/max(max(plots(:,:,2))) );
    
    solRank = find(bestRuns(:,1) == ii);
    
    
    I = imtile(flipud(plots));
    imagesc(I)
    axis off
    caxis([minVal maxVal]);
    
    title(['Run ', num2str(ii),', Rank ',num2str(solRank)])
end

%% All Sols Diffs Plot

clear plots I 
figure(98)

p = numSubPlots(data.numRuns);

for ii = 1:data.numRuns
    subplot(p(1),p(2),ii)
    
    plots(:,:,1) = abs(beam_orig.field_fList).^2 -...
        abs(theHerd(ii).field_fList).^2;
    plots(:,:,2) = abs( angle(beam_orig.field_fList) - ...
        angle(theHerd(ii).field_fList) );
    plots(:,:,2) = max(max(plots(:,:,1))) * (plots(:,:,2)/max(max(plots(:,:,2))) );  
    
    solRank = find(bestRuns(:,1) == ii);
    
    I = imtile(flipud(plots));
    imagesc(I)
    axis off
    colorbar
%     caxis([minVal maxVal]);
    
    title(['Run ', num2str(ii),', Rank ',num2str(solRank)])
end

%% resave the data (seriously a bad idea lol)

a = fieldnames(data);

for ii = 1:length(a)
    
    eval([a{ii},' = data.',a{ii},';']);
    
end

save('C:\Users\rlemons\Documents\finalBeamsNEW.mat',a{:})

%% Testbed

%% Suport Functions
function [beam_orig,beamsHerd] = reconBeams(data)

beam_orig = beamPropagation2D(data.beam_orig_BAK);
[beam_orig.field_fList,beam_orig.field_FList] = ...
    beam_orig.forwardProp_FreeSpace2D(data.z);

herdParams = beam_orig.outputProperties2D(data.beamType);


for ii = 1:data.numRuns
    beamsHerd(ii) = beamPropagation2D(data.beam_orig_BAK);
end

for ii = 1:data.numRuns
    for jj = 1:size(data.propTest,1)
        herdParams.([data.beamType,'_',data.propTest{jj,1}]) = ...
            data.sol(ii).finalProps.([data.beamType,'_',data.propTest{jj,1}])(end,:);
    end
    beamsHerd(ii).([data.beamType,'_','InitialBeamDef2D'])(herdParams);
    [beamsHerd(ii).field_fList,beamsHerd(ii).field_FList] = ...
        beamsHerd(ii).forwardProp_FreeSpace2D(data.z);
end

if ~strcmpi(data.solveType,'beam')
    beam_orig = double(imread(data.fileName));
    beam_orig = beam_orig/max(max(beam_orig));
end

end

function p = numSubPlots(n)

p(2) = ceil(sqrt(n));
p(1) = ceil(n/p(2));

end
