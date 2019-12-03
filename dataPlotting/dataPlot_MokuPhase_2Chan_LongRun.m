% General file to plot a number of different possible aspects. Comment and
% uncomment as needed.

%% Run this section find the files
clear;
% close all;

% Read in files
[fileName,filePath] = uigetfile('D:\CEP\CSVs\*.*');
fileStr = fullfile(filePath,fileName);

% Turn off warning about number of peaks. QOL change, can be removed.
warning('off','signal:findpeaks:largeMinPeakHeight');

%% Load the data

% Changed up so the data is only loaded in to memory once
dataOUT = dataImport_Moku(fileStr);
numChans = (size(dataOUT,2)-1)/5;

% With this indexing you get frequency, phase, I, Q, amplitude
% Usage = vars(chan,variable) ie. vars(2,3) = row 2, column 3 = Chan 2's I
chans = {'chan1','chan2'};
vars = [2 3 4 5 6;7 8 9 10 11];


%%%%% Scale Time %%%%%
dataOUT(:,1) = dataOUT(:,1)/60/60; % Choose unit of time ie. seconds/60/60 = hours


%%%%% Remove Set frequency %%%%%
for ii = 1:length(chans)
    setFreq.(chans{ii}) = dataOUT(1,vars(ii,1));
end
dataOUT(:,2:5) = dataOUT(:,3:6); % Shift Columns

if length(chans) == 2
    dataOUT(:,7:10) = dataOUT(:,8:11);
end

%%%%% Remove huge frequency outliers %%%%%
for ii = 1:length(chans)
    out = abs((dataOUT(:,vars(ii,1))-mean(dataOUT(:,vars(ii,1)))))...
        > 5*std(dataOUT(:,vars(ii,1)));
    dataOUT(out,vars(ii,1)) = mean(dataOUT(:,vars(ii,1)));
end

%%%%% Unwrap large phase jumps %%%%%
for ii = 1:length(chans)
    dataOUT(:,vars(ii,2)) = unwrap(dataOUT(:,vars(ii,2)),10^6);
end


clear out

%% Process the data


% Clear and pre-allocate phase array info
clear vals phaseOUT dateAndTimes
vals = zeros(length(chans),2);
phaseOUT = zeros(length(dataOUT(:,1)),length(chans));

startTime = datetime(2019,08,15,11,20,00);
endTime = startTime + hours(24);
dateAndTimes = linspace(startTime,endTime,length(dataOUT(:,1)))';

%%%%% Perform Operations for Each Channel %%%%%
for ii = 1:numChans % second loop to gather the second channel data
    
    %%%%% Frequency Section %%%%%
    
    % Average Frequency
    meanVal.(chans{ii}).freq = mean(dataOUT(:,vars(ii,1)));
    minVal.(chans{ii}).freq = min(dataOUT(:,vars(ii,1)));
    maxVal.(chans{ii}).freq = max(dataOUT(:,vars(ii,1)));
    span.(chans{ii}).freq =...
        maxVal.(chans{ii}).freq - minVal.(chans{ii}).freq;
    
    
    %%%%% Phase Section %%%%%
    
    % Average Phase
    meanVal.(chans{ii}).phase = mean(dataOUT(:,vars(ii,2)));
    minVal.(chans{ii}).phase = mean(dataOUT(:,vars(ii,2)));
    maxVal.(chans{ii}).phase = mean(dataOUT(:,vars(ii,2)));
    span.(chans{ii}).phase =...
        maxVal.(chans{ii}).phase - minVal.(chans{ii}).phase;
    
    % Fit to phase info to remove drifts due to offset in Moku
    vals(ii,:) = polyfit(dataOUT(:,1),dataOUT(:,vars(ii,2)),1);
    phaseOUT(:,ii) = dataOUT(:,vars(ii,2))...
        - vals(ii,1).*dataOUT(:,1) + vals(ii,2);
    
    phaseOUT(:,ii) = phaseOUT(:,ii) - mean(phaseOUT(:,ii)); % Center phase info around 0
    phaseOUT(:,ii) = phaseOUT(:,ii) .* 6.2832; % Convert from cycles to radians
    
    
    %%%%% Amplitude Section %%%%%
    
    dataOUT(:,vars(ii,5)) =...
        sqrt(dataOUT(:,vars(ii,3)).^2+dataOUT(:,vars(ii,4)).^2); % add I and Q in quadrature to get A
    
    % Average Amplitude
    meanVal.(chans{ii}).amp = mean(dataOUT(:,vars(ii,5)));
    minVal.(chans{ii}).amp = min(dataOUT(:,vars(ii,5)));
    maxVal.(chans{ii}).amp = max(dataOUT(:,vars(ii,5)));
    span.(chans{ii}).amp =...
        maxVal.(chans{ii}).amp - minVal.(chans{ii}).amp;
    
end

%% Allan Varience of time signal

for jj = 1:numChans
    
    alStruc(jj).freq = dataOUT(:,vars(jj,1))/meanVal.(chans{jj}).freq-1;
    alStruc(jj).rate = length(dataOUT(:,1))/24/60/60
    
    tmp = (1/alStruc(jj).rate)*logspace(0,6,50);
    tmp(tmp > 10^4) = [];
    
    alStruc(jj).taus = tmp;
    
end

alStruc(1).name = 'f_{CEO}';
alStruc(2).name = 'f_{LO}';

for jj = 1:2
    allan(alStruc(jj),alStruc(jj).taus,alStruc(jj).name);
end

%% Phase Noise

% pVec = zeros(size(dataOUT(:,1),2),2);
% 
% for ii = 1:numChans
%     pVec(:,ii) = cumtrapz(dataOUT(:,1) , (dataOUT(:,vars(ii,1))-meanVal.(chans{ii}).freq));
% end

N = length(dataOUT(:,1));

sampRate = length(dataOUT(:,1))/24/60/60;
freqVec = sampRate*((0:N/2)/N);

ii = 2;
fftTime = fft( dataOUT(:,vars(ii,1))- meanVal.(chans{ii}).freq );
fftTime = abs(fftTime/N);
fftTime = fftTime(1:N/2+1);
fftTime(2:end-1) = 2*fftTime(2:end-1);

fftPhase = fftTime' .*( (meanVal.(chans{ii}).freq^2)./freqVec.^2);

hold on
loglog(freqVec,fftTime)
hold off


%% Run this section to plot the data

chansPlot = 1:numChans;
avg = 2;
numTicks = 8;

%%% Chose Which parameters to plot %%%
% [frequency, phase, amplitude] %
plotParams = [1, 0, 0];


%%%%% Plot Frequency %%%%%
if plotParams(1) == 1
    
    
    titles = {'Frequency Drift $f_{AOFS}$, w/ longterm PID control',...
        sprintf('Frequency Drift, Averaging = %d, $f_{OOL}$',avg)};
    %     xLabel = repmat({'Time (h)'},1,numChans);
    xLabel = repmat({'Time (h)'},1,numChans);
%     yLabel = repmat({'$$\Delta \nu$$ (Hz)'},1,numChans);
    yLabel = {'$$\Delta \nu$$ (kHz)','$$\Delta \nu$$ (Hz)'};
    
    yScale = [1e3,1e0];
    lineW = [1,1];
    
    
    
    for ii = chansPlot
        figure(ii);
        ax = gca;
        
        plot(dataOUT(:,1),...
            movmean(...
            (dataOUT(:,vars(ii,1))-meanVal.(chans{ii}).freq)/yScale(ii),...
            avg),...
            'LineWidth', lineW(ii),'color',[33 54 86]/255);

%         plot(dateAndTimes,...
%             movmean(...
%             (dataOUT(:,vars(ii,1))-meanVal.(chans{ii}).freq)/yScale(ii),...
%             avg),...
%             'LineWidth', lineW(ii),'color',[33 54 86]/255);        
        
        rmsFreq = rms(...
            movmean(...
            (dataOUT(:,vars(ii,1))-meanVal.(chans{ii}).freq)/yScale(ii),...
            avg)...
            );
        
        disp(['rms Freq Jitter: ',num2str(rmsFreq)])
        
        
        

%         xlim([dateAndTimes(1) dateAndTimes(end)]);
%         xtickangle(60);
%         xticks(round(dataOUT(1:floor(length(dataOUT(:,1))/numTicks):end,1)));
        xlim([min(dataOUT(:,1)) max(dataOUT(:,1))]);

        
        
        % Make it look nice
        xlabel(xLabel{ii},'FontSize',24,'Interpreter','latex');
        ylabel(yLabel{ii},'FontSize',24,'Interpreter','latex');
        ax.FontSize = 60;
%         title(titles{ii},'FontSize',40,'Interpreter','latex');
        %         legend({'No Averaging','Averaging = 5'},'Location','northeast');
    end
    
    
end


%%%%% Plot Phase %%%%%
if plotParams(2) == 1
    
    
    
    titles = {sprintf('Phase Drift, Averaging = %d, $f_{AOFS}$',avg),...
        sprintf('Phase Drift, Averaging = %d, $f_{OOL}$',avg)};
    %     xLabel = repmat({'Time (h)'},1,numChans);
    xLabel = repmat({'Time (min)'},1,numChans);
    yLabel = repmat({'Phase (rad)'},1,numChans);
    lineW = [1,1];
    
    
    
    for ii = chansPlot
        
        figure(ii+numChans);
        ax = gca;
        
        plot(dataOUT(:,1),...
            movmean(phaseOUT(:,ii),avg),...
            'LineWidth', lineW(ii));
        
        
        rmsPhase = rms(...
            movmean(phaseOUT(:,ii),avg)...
            );
        
        disp(['rms Phase Jitter: ',num2str(rmsPhase)])
        
        
        
        xlim([min(dataOUT(:,1)) max(dataOUT(:,1))]);
        
        
        % Make it look nice
        xlabel(xLabel{ii},'FontSize',24,'Interpreter','latex');
        ylabel(yLabel{ii},'FontSize',24,'Interpreter','latex');
        ax.FontSize = 34;
        title(titles{ii},'FontSize',40,'Interpreter','latex');
        %         legend({'No Averaging','Averaging = 5'},'Location','northeast');
        
    end
    
end


%%%%% Plot Amplitude %%%%%
if plotParams(3) == 1
    
    titles = {sprintf('Amplitude Drift, Averaging = %d, $f_{AOFS}$',avg),...
        sprintf('Amplitude Drift, Averaging = %d, $f_{OOL}$',avg)};
    %     xLabel = repmat({'Time (h)'},1,numChans);
    xLabel = repmat({'Time (min)'},1,numChans);
    yLabel = repmat({'Amplitude (mV)'},1,numChans);
    lineW = [1,1];
    
    yScale(1:2) = 10^-3;
    
    
    for ii = chansPlot
        figure(ii+numChans*2);
        ax = gca;
        
        plot(dataOUT(:,1),...
            movmean(...
            (dataOUT(:,vars(ii,5)))/yScale(ii),...
            avg),...
            'LineWidth', lineW(ii),'color',[33 54 86]/255);
        
        
        rmsFreq = rms(...
            movmean(...
            (dataOUT(:,vars(ii,5))-meanVal.(chans{ii}).amp)/yScale(ii),...
            avg)...
            );
        
        disp(['Amplitude Jitter: ',num2str(rmsFreq)])
        
        
        
        
        xlim([min(dataOUT(:,1)) max(dataOUT(:,1))]);
        
        
        % Make it look nice
        xlabel(xLabel{ii},'FontSize',24,'Interpreter','latex');
        ylabel(yLabel{ii},'FontSize',24,'Interpreter','latex');
        ax.FontSize = 40;
        title(titles{ii},'FontSize',40,'Interpreter','latex');
        %         legend({'No Averaging','Averaging = 5'},'Location','northeast');
    end
    
    
end

%% Special Plots %%

avg = 1;

%%%%% Phase Subtraction Plot %%%%%

titles = {'Frequency Drift, LOCKED minus UNLOCKED'};
xLabel = repmat({'Time (h)'},1,numChans);
yLabel = repmat({'Phase (rad)'},1,numChans);

figure(99);
% subplot(3,1,3);
ax = gca;

plot(dataOUT(:,1),...
    movmean(phaseOUT(:,1)-phaseOUT(:,2),avg)...
    );
rmsPhase = rms(...
    movmean(phaseOUT(:,1)-phaseOUT(:,2),avg)...
    );
disp(['rms Freq Jitter: ',num2str(rmsPhase)])


xlim([min(dataOUT(:,1,1)) max(dataOUT(:,1,1))]);
% legend({'LOCKED','UNLOCKED'},'Location','northeast');


% Make it look nice
xlabel(xLabel{1},'FontSize',24,'interpreter','latex');
ylabel(yLabel{1},'FontSize',24,'interpreter','latex');
ax.FontSize = 34;
title(titles{1},'FontSize',40,'interpreter','latex');


%% Special Data Analysis
% Stuff that really only needs to be done a couple times. Not important or
% doesn't make enough sense to use in main loops


%%% Special time range data analysis of Phase data
lowTime = 5;
[~,shiftLow] = min(abs(dataOUT(:,1,1)-lowTime));
highTime = 5.2;
[~,shiftHigh] = min(abs(dataOUT(:,1,1)-highTime));

if shiftLow ~= 1 && shiftHigh ~= length(dataOUT(:,1,1))
    phaseNew = phaseOUT(shiftLow:shiftHigh,:);
end


for ii = 1:size(phaseNew,2)
    
    phaseNew(:,ii) = phaseNew(:,ii) - mean(phaseNew(:,ii));
    
    figure(50-ii)
    avg = 10;
    plot(dataOUT(shiftLow:shiftHigh,1),movmean(phaseNew(:,ii),avg))
    xlim([dataOUT(shiftLow,1) dataOUT(shiftHigh,1)]);
    disp([sprintf('rms phase during %.1fh - %.1fh: ',lowTime,highTime),...
        num2str(rms(phaseNew(:,ii)))])
end