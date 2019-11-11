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
nm = 10^-9 * m;
km = 10^3 * m;

z =  40 * m; % Distance to Propagate
lambda = 1.5 * um; % Wavelength of light

plotFlag = 1;

numHerd = 30;
numGens = 60;
numRuns = 1;

% tst = mfilename('fullpath');
% tst = dbstack;

beamType = 'hex';
propTest = {'PhaseOffset',[0,2*pi],0.1};
% propTest = {'PhaseOffset',[0,2*pi];'AmpBeams',[0,1]};
% propTest = {'PhaseOffset',[0,2*pi],0.1;...
%     'AmpBeams',[0,1],0.1;...
%     'BeamsOn',[0,1],0.1;};

for ii = 1:numel(propTest(:,1))
    propList{ii,1} = [beamType, '_', propTest{ii,1}];
end

finalFileName = 'C:\Users\rlemons\Documents\GitHub\varCodes\savedVars\finalBeams.mat';
beamFileName = 'C:\Users\rlemons\Documents\GitHub\varCodes\savedVars\genAlg_Beams.mat';


%% Generate and save Patient Zero
beam_orig = beamPropagation2D(lambda,6*cm,6*cm,2^10,'hex');

save(beamFileName,'beam_orig');

%% Load Patient Zero
load(beamFileName);
beam_orig_BAK = beam_orig.outputProperties2D('');
herdParams = repmat(beam_orig.outputProperties2D(beamType),numHerd,1);

% fileName = 'C:\Users\rlemons\Documents\tst.png';
fileName = [];

if isempty(fileName)
% Prop the original beam to create the problem to solve
[beam_orig.field_fList,beam_orig.field_FList] =...
    beam_orig.forwardProp_FreeSpace2D(z);
minVal = min(min(abs(beam_orig.field_fList).^2));
maxVal = max(max(abs(beam_orig.field_fList).^2));

solveType = 'beam';

save(finalFileName,'solveType','beam_orig_BAK','z',...
    'numHerd','numGens','numRuns','beamType','propTest');

elseif contains(fileName,'.mat')
    
    error('Haven''t coded this yet... My B');
    
%     solveType = 'matrix';
%     
%     save(finalFileName,'solveType','fileName''beam_orig_BAK','z',...
%         'numHerd','numGens','numRuns','beamType','propTest');
    
else
    
    [solImage,cmap,~] = imread(fileName);
    
    if ~isempty(cmap)
        solImage = rgb2gray(ind2rgb(solImage,cmap));
    elseif size(solImage,3) == 3
        solImage = rgb2gray(solImage);
    end
    
    solImage = double(...
        imresize(solImage,[beam_orig.grid_npts beam_orig.grid_npts])...
        );
    solImage = solImage./max(max(solImage));
    
    minVal = min(min(solImage));
    maxVal = max(max(solImage));
    
    solveType = 'image';
    
    save(finalFileName,'solveType','fileName','beam_orig_BAK','z',...
        'numHerd','numGens','numRuns','beamType','propTest');
    
end

for zz = 1:numRuns
    
    %%%%% The Herd
    
    clear theHerd herdProps
    
    
    zTime = tic;
    for ii = 1:numHerd
        theHerd(ii) = beamPropagation2D(beam_orig_BAK);
    end
    
    for jj = 1:numel(propList)
        
        herdProps.(propList{jj}) = repmat(beam_orig.(propList{jj}),numHerd,1);
        
        for ii = 1:numHerd
            
            herdProps.(propList{jj})(ii,:) =...
                createInd(propTest{jj,2},numel(herdProps.(propList{jj})(ii,:)),propList{jj});
            
        end
        
    end
    
    herdParams = importProps(herdParams,herdProps);
    
    
    for ii = 1:numHerd
        
        theHerd(ii).([beamType,'_InitialBeamDef2D'])(herdParams(ii));
        
    end
    
    errAVG = zeros(1,numGens);
    for ii = 1:numel(propList)
        propsFinal.(propList{ii}) = zeros(numGens,numel(herdProps.(propList{ii})(1,:)));
    end
    
    if plotFlag
        subplot(1,2,1); %#ok<*UNRCH>
        if strcmpi(solveType,'beam')
            beam_orig.plotField2D(beam_orig.field_fList,'abs');
            caxis([minVal maxVal]);
        elseif strcmpi(solveType,'matrix')
        else
            imagesc(solImage);
            axis square off
%             caxis([minVal maxVal]);
        end
        subplot(1,2,2);
    end
    
    for ii = 1:numGens
        
        gTime = tic;
        
        for jj = 1:numHerd
            %%%%% Propagate
            [theHerd(jj).field_fList,theHerd(jj).field_FList] =...
                theHerd(jj).forwardProp_FreeSpace2D(z);
        end
        
        %%%%% Evaluate: The Herd
        if strcmpi(solveType,'beam')
            [err,errAVG(ii)] = judge(beam_orig,theHerd); % Report Error
        elseif strcmpi(solveType,'matrix')
            
        else
            [err,errAVG(ii)] = judge(solImage,theHerd); % Report Error
        end
        theHerd = theHerd(err(:,1)); % Sort by best (Decending)
        
        sol(zz).fitness(ii,1) = err(1,2);
        sol(zz).fitness(ii,2) = errAVG(ii);
        
        for jj = 1:numel(propList)
            propsFinal.(propList{jj})(ii,:) = theHerd(1).(propList{jj});
        end
        
        disp(['>> Error of Gen ',num2str(ii),' is: ',num2str(errAVG(ii),'%.1f\n')])
        disp(['>> The best individual is: ',num2str(err(1,1))])
        
        if plotFlag
            theHerd(1).plotField2D(theHerd(1).field_fList,'abs')
%             caxis([minVal maxVal]);
            drawnow
            pause(0.5)
        end
        
        if ii < numGens
            %%%%% Recreate: The Herd
            theHerd = cull(theHerd,beamType,propTest);
        end
        sol(zz).timeG(ii) = toc(gTime);
    end
    
    sol(zz).timeZ = toc(zTime);    
%     sol(zz).beam = theHerd(1);
    sol(zz).finalProps = propsFinal;
    save(finalFileName,'sol','-append')
    
end

fclose('all');

%% Functions for use in the calc

function genes = createInd(range,nBeams,str)

if contains(str,'PhaseOffset')
    genes = rand(nBeams,1).*(range(2) - range(1)) + range(1);
    genes(1) = 0;
elseif contains(str,'AmpBeams')
    genes = rand(nBeams,1).*(range(2) - range(1)) + range(1);
elseif contains(str,'BeamsOn')
    genes = randi(2,nBeams,1)-1;
elseif contains(str,'PhaseCurve')
    genes = repmat(rand.*(range(2) - range(1)) + range(1),nBeams,1);
end


end

function params = importProps(params,props)

list = fieldnames(props);

for ii = 1:numel(list)
    for jj = 1:numel(params)
        params(jj).(list{ii}) = props.(list{ii})(jj,:);
    end
end

end

function [err,errAvg] = judge(ideal,herd)

numHerd = numel(herd);

err = zeros(numHerd,2);
err(:,1) = 1:numHerd;
for ii = 1:numHerd
    
    %%%%% Error generation
    if isa(ideal,'beamPropagation2D')
%         err(ii,2) = numel(ideal.field_fList)/...
%             abs(sum(abs(ideal.field_fList-herd(ii).field_fList).^2,'all'));
        err(ii,2) = numel(ideal.field_fList)/...
            abs(sum(abs( ...
            ideal.field_fList./max(max(ideal.field_fList)) -...
            herd(ii).field_fList./max(max(ideal.field_fList)) ...
            ).^2,'all'));
    else
%         err(ii,2) = numel(ideal)/abs(sum((ideal-abs(herd(ii).field_fList).^2).^2,'all'));
        err(ii,2) = 1/(1-slimSSIM(ideal,abs(herd(ii).field_fList).^2));
    end
    
%     if err(ii,2) == Inf
%         disp('>> The values for an exact solution are:')
%         fprintf('\n\n')
%         disp(num2str(herd(1).hex_PhaseOffset,'%0.4f'))
%         fprintf('\n\n')
%         error('You are good son!');
%     end
    
end
err = sortrows(err,-2);
errAvg = mean(err(:,2));

end

function herd = cull(herd,beamType,propList)

N = numel(herd);
numBest = floor(N * .14);
numRand = ceil(N/5) - numBest;
numKeep = numBest+numRand;

randKeep = randperm(length(numBest+1:N),numRand)+numBest;
herd(numBest+1:numKeep) = herd(randKeep);


count = 0;
while count < N - numKeep
    mom = randi(numKeep);
    dad = randi(numKeep);
    if mom ~= dad
        herd(numKeep+count+1) =...
            breedParents(herd(mom),herd(dad),beamType,propList);
        count = count + 1;
    end
end

for ii = 1:N
    herd(ii).([beamType,'_InitialBeamDef2D'])(herd(ii).outputProperties2D(beamType));
end

end

function baby = breedParents(mom,dad,beamType,propList)

if rand > 0.5
    nCross = ceil(dad.([beamType,'_NBeams'])/2);
else
    nCross = floor(dad.([beamType,'_NBeams'])/2);
end

dadGenes = randperm(dad.([beamType,'_NBeams'])-1,nCross)+1;

baby = beamPropagation2D(mom);

for ii = 1:size(propList,1)
    
    propStr = [beamType,'_',propList{ii,1}];
    
    baby.(propStr)(dadGenes) = dad.(propStr)(dadGenes);
    
    if contains(propStr,'PhaseOffset')
        for jj = 2:length(baby.(propStr))
            if propList{ii,3} > rand
                baby.(propStr)(randi(baby.([beamType,'_NBeams'])-1)+1) = ...
                    rand.*(propList{ii,2}(2)-propList{ii,2}(1)) + propList{ii,2}(1);
            end
        end
    elseif contains(propStr,'AmpBeams')
        for jj = 1:length(baby.(propStr))
            if propList{ii,3} > rand
                baby.(propStr)(randi(baby.([beamType,'_NBeams']))) = ...
                    rand.*(propList{ii,2}(2)-propList{ii,2}(1)) + propList{ii,2}(1);
            end
        end
    elseif contains(propStr,'BeamsOn')
        for jj = 1:length(baby.(propStr))
            if propList{ii,3} > rand
                baby.(propStr)(randi(baby.([beamType,'_NBeams']))) = ...
                    randi(2)-1;
            end
        end
    elseif contains(propStr,'PhaseCurve')
        if propList{ii,3} > rand
            baby.(propStr) = repmat(...
                rand.*(propList{ii,2}(2)-propList{ii,2}(1)) + propList{ii,2}(1),...
                baby.([beamType,'_NBeams']),1);
        end
    end
    
end

end








