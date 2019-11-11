%% Main file for plotting data from Optical Spec Analyzer.
% Comment and uncomment as needed

%% Load the data
clear;
close all;

% Read in files
[dataIN,fileStr] = dataImport_OOS(); %#ok<*SAGROW>
numFilesIn = size(dataIN,3);



%% Process the data

dataOUT = dataIN;

%%% Select ranges for data analysis
lowPass = 520; % in nm
highPass = 530; % in nm

[~,shift] = min(abs(dataOUT(:,1,1)-lowPass));
if shift ~= 1
    dataOUT(1:shift,:,:) = [];
end
[~,shift] = min(abs(dataOUT(:,1,1)-highPass));
if shift ~= length(dataOUT(:,1,1))
    dataOUT(shift:end,:,:) = [];
end

%%% Interpolate data to twice the spaceing
expandFac = 3;
tmp = zeros(numel(dataOUT(:,1,1))*expandFac,2,numFilesIn);
for ii = 1:numFilesIn
    tmp(:,1,ii) = linspace(...
        dataOUT(1,1,ii),dataOUT(end,1,ii),...
        numel(dataOUT(:,1,ii))*expandFac...
        );
    tmp(:,2,ii) = interp1(...
        dataOUT(:,1,ii),dataOUT(:,2,ii),tmp(:,1,ii)...
        );
end
dataOUT = tmp;


for ii = 1:numFilesIn
    
    % Important Values
    minVal(ii).counts = min(dataOUT(:,2,ii));
    [maxVal(ii).counts,pos(ii)] = max(dataOUT(:,2,ii));
    
    % Normalize the data
    dataOUT(:,2,ii) = dataOUT(:,2,ii)/maxVal(ii).counts;
    
    % FWHM
    ind1 = find(dataOUT(:,2,ii) >= 0.5,1,'first');
    ind2 = find(dataOUT(:,2,ii) >= 0.5,1,'last');
    fwhm(ii,1) = dataOUT(ind2,1,ii) - dataOUT(ind1,1,ii);
    
end


%% Plot the data

numPlot = 1:size(dataOUT,3);

%%%% Individual File Plotting %%%%
xLabels = repmat({'Wavelength (nm)'},numFilesIn,1);
yLabels = repmat({'Intensity (arb.)'},numFilesIn,1);


for ii = numPlot
    figure(ii);
    ax = gca;
    
    [~,titles,~] = fileparts(fileStr{ii});
    
    
    plot(dataOUT(:,1,ii),dataOUT(:,2,ii),...
        'linewidth',3,'color',[33 54 86]/255);
    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    
    
    
    title(titles,'FontSize',30,'Interpreter','none');
    xlabel(xLabels{ii},'FontSize',24);
    ylabel(yLabels{ii},'FontSize',24);
    
    ax.FontSize = 44;
    
end

%% Special Plots

numPlot = 1:size(dataOUT,3);

%%%% All plot on one %%%%
titles = 'All Spectra Plot';
xLabels = 'Wavelength (nm)';
yLabels = 'Intensity (arb.)';

figure(99);
ax = gca;

legendVals = split(fileStr,'_');
legendVals = char(erase(legendVals(:,:,end)','.txt'));
legendVals = [legendVals,...
    repmat(', FWHM: ',numFilesIn,1),...
    num2str(fwhm,'%-.2f')];


hold off
for ii = numPlot
    plot(dataOUT(:,1,ii),dataOUT(:,2,ii),...
        'linewidth',1);
    hold on
end
hold off

xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);

legend(legendVals,'Location','northeast');

title(titles,'FontSize',30,'Interpreter','none');
xlabel(xLabels,'FontSize',24);
ylabel(yLabels,'FontSize',24);

ax.FontSize = 44;



%% Empty space because i hate looking at the bottom of my screen.
% I mean, really, who wants to have to look at multiple places on their
% screen instead of just being able to continue scrolling down even if
% those lines are empty. That is quite possibly the stupidest way of
% setting up your IDE that I am baffeled. They all seem to do it and I KNOW
% I can't be the only person who just wants to scroll down.
%
%
%
% Really I only need 10 extra lines you see. Nothing crazy. But NO, I'm the
% crazy one apparently.