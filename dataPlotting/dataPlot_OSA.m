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
    
    minVal(ii).dBm = min(dataOUT(:,2,ii));
    [maxVal(ii).dBm,pos(ii)] = max(dataOUT(:,2,ii)); 
    
end




%% Plot the data

numPlot = 1:size(dataOUT,3);


titles = {'In Loop Sprectrum',...
    'Out Of Loop Spectrum'};
xLabels = repmat({'Wavelength (nm)'},numFilesIn+1,1);
yLabels = repmat({'Power (dBm)'},numFilesIn+1,1);



for ii = numPlot
    figure(ii);
    ax = gca;
    
    
    plot(dataOUT(:,1,ii),dataOUT(:,2,ii));
    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    ylim([minVal(ii).dBm-10 maxVal(ii).dBm+10]);
    
    ax.FontSize = 24;
    
    title(titles{ii},'FontSize',30);
    xlabel(xLabels{ii},'FontSize',24);
    ylabel(yLabels{ii},'FontSize',24);
    
    
    
end