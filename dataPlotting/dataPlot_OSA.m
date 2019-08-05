%% Main file for plotting data from Optical Spec Analyzer.
% Comment and uncomment as needed

%% Load the data
clear;
close all;

% Read in files
[dataIN,fileStr] = dataImport_OSA(); %#ok<*SAGROW>
numFilesIn = size(dataIN,3);



%% Process the data

dataOUT = dataIN;


for ii = 1:numFilesIn

    % Average the data
    avgVal = 5;
    dataOUT(:,2,ii) = movmean(dataOUT(:,2,ii),avgVal);
    
    % Force Ceil or Floor for data
    floorForce = -75;
    dataOUT( dataOUT(:,2,ii)<floorForce ,2,ii) = floorForce;
    ceilForce = 0;
    dataOUT( dataOUT(:,2,ii)>ceilForce ,2,ii) = ceilForce;
        
    % Important Values
    minVal(ii).dBm = min(dataOUT(:,2,ii));
    [maxVal(ii).dBm,pos(ii)] = max(dataOUT(:,2,ii));
    
end




%% Plot the data

numPlot = 1:size(dataOUT,3);


titles = {'Out Of Loop Sprectrum',...
    'Out Of Loop Spectrum'};
xLabels = repmat({'Wavelength (nm)'},numFilesIn+1,1);
yLabels = repmat({'Power (dBm)'},numFilesIn+1,1);



for ii = numPlot
    figure(ii);
    ax = gca;
    
    
    plot(dataOUT(:,1,ii),dataOUT(:,2,ii),...
        'linewidth',3,'color',[33 54 86]/255);
    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    ylim([minVal(ii).dBm maxVal(ii).dBm+10]);
    

    
    title(titles{ii},'FontSize',30);
    xlabel(xLabels{ii},'FontSize',24);
    ylabel(yLabels{ii},'FontSize',24);
    
    ax.FontSize = 44;
        
end