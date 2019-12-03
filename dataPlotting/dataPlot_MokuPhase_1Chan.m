% General file to plot a number of different possible aspects. Comment and
% uncomment as needed.

%% Run this section to load the data
clear;
close all;

% Read in files
[dataIN,fileStr] = dataImport_Moku(); %#ok<*SAGROW>
numFilesIn = size(dataIN,3);

% Turn off warning about number of peaks. QOL change, can be removed.
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Run this section to process the data and find peaks

dataOUT = dataIN;

% % Offset for any attenuation
% dataOUT(:,2,:) = dataOUT(:,2,:)+10;
%
% %%% Select ranges for data analysis
% highPass = 31.1e6;
% [~,shift] = min(abs(dataOUT(:,1,1)-highPass));
% if shift ~= 1
%     dataOUT(1:shift,:,:) = [];
% end
%
% lowPass = 31.8e6;
% [~,shift] = min(abs(dataOUT(:,1,1)-lowPass));
% if shift ~= length(dataOUT(:,1,1))
%     dataOUT(shift:end,:,:) = [];
% end



%%% Perform operations on each file
for ii = 1:numFilesIn
    
    %%%%% Time Section %%%%%
    
    dataOUT(:,1,ii) = dataOUT(:,1,ii)/60/60;
    
    %%%%% Frequency Section %%%%%
    
    % Average Frequency
    meanVal.freq(ii) = mean(dataOUT(:,3,ii));
    minVal.freq(ii) = min(dataOUT(:,3,ii));
    maxVal.freq(ii) = max(dataOUT(:,3,ii));
    span.freq(ii) = maxVal.freq(ii)-minVal.freq(ii);
    
    
    %%%%% Phase Section %%%%%
    
    % Average Phase, (meaningless for one signal)
    meanVal.phase(ii) = mean(dataOUT(:,4,ii));
    minVal.phase(ii) = min(dataOUT(:,4,ii));
    maxVal.phase(ii) = max(dataOUT(:,4,ii));
    span.phase(ii) = maxVal.phase(ii)-minVal.phase(ii);
    
    %%%%% Amplitude Section %%%%%
    
    dataOUT(:,7,ii) = sqrt(dataOUT(:,5,ii).^2+dataOUT(:,6,ii).^2);
    
    % Average Amplitude
    meanVal.amp(ii) = mean(dataOUT(:,7,ii));
    minVal.amp(ii) = min(dataOUT(:,7,ii));
    maxVal.amp(ii) = max(dataOUT(:,7,ii));
    span.amp(ii) = maxVal.amp(ii)-minVal.amp(ii);
    
end



%% Run this section to plot the data

numPlots = numFilesIn;
% numPlots = size(dataOUT,2);
avg = 5;


%%%%% Plot Frequency %%%%%

titles = {'Frequency Drift, LOCKED',...
    'Frequency Drift, UNLOCKED'};
xLabels = repmat({'Time (h)'},numPlots,1);
% yLabels = repmat({'Frequency (MHz)'},numPlots,1);
yLabels = {'$$\Delta \nu$$ (kHz)',...
    '$$\Delta \nu$$ (kHz)'};
yScale = [1e3,1e3];

% figure(1);

for ii = 1:numPlots

    figure(ii);
%     subplot(3,1,ii)
    ax = gca;
    
    plot(dataOUT(:,1,ii),...
        movmean((dataOUT(:,3,ii)-meanVal.freq(ii))/yScale(ii),avg),...
        'LineWidth',2);
    
    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    ylim([...
        ((minVal.freq(ii)-(span.freq(ii)*.25/2))-meanVal.freq(ii))/yScale(ii)...
        ((maxVal.freq(ii)+(span.freq(ii)*.25/2))-meanVal.freq(ii))/yScale(ii)...
        ]);
    ax.FontSize = 24;
    
    % Make it look nice
    title(titles{ii},'FontSize',30,'Interpreter','latex');
    xlabel(xLabels{ii},'FontSize',24,'Interpreter','latex');
    ylabel(yLabels{ii},'FontSize',24,'Interpreter','latex');
    
    
    
end

%%%%% Plot Phase %%%%%

titles = {'Long Term Phase Drift, LOCKED',...
    'Long Term Phase Drift, UNLOCKED'};
xLabels = repmat({'Time (min)'},numPlots,1);
yLabels = repmat({'!!!!ADD THIS!!!!'},numPlots,1);


for jj = 1:numPlots

    figure(jj+numFilesIn);
    ax = gca;

    plot(dataOUT(:,1,jj),dataOUT(:,3,jj));

    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    ylim([...
        (minVal.freq(ii) - (span.freq(ii)*.2/2))/1e6...
        (maxVal.freq(ii) + (span.freq(ii)*.2/2))/1e6...
        ]);
    ax.FontSize = 24;

    % Make it look nice
    title(titles{jj},'FontSize',40,'interpreter','latex');
    xlabel(xLabels{jj},'FontSize',34,'interpreter','latex');
    ylabel(yLabels{jj},'FontSize',34,'interpreter','latex');

end

%%%%% Plot Amplitude %%%%%

titles = {'Long Term Amplitude Drift, LOCKED',...
    'Long Term Amplitude Drift, UNLOCKED'};
xLabels = repmat({'Time (min)'},numPlots,1);
yLabels = repmat({'Amplitude (V)'},numPlots,1);


for kk = 1:numPlots

    figure(kk+(2*numFilesIn));
    ax = gca;

    plot(dataOUT(:,1,kk),dataOUT(:,3,kk));

    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    ylim([...
        (minVal.freq(ii) - (span.freq(ii)*.2/2))/1e6...
        (maxVal.freq(ii) + (span.freq(ii)*.2/2))/1e6...
        ]);
    ax.FontSize = 24;

    % Make it look nice
    title(titles{kk},'FontSize',40,'interpreter','latex');
    xlabel(xLabels{kk},'FontSize',34,'interpreter','latex');
    ylabel(yLabels{kk},'FontSize',34,'interpreter','latex');

end

%% Special Plots %%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(99);
% 
% titles = {'Frequency Drift, $$f_{OOL}$$ and $$f_{CEP}$$'};
% xLabels = repmat({'Time (min)'},numPlots,1);
% yLabels = repmat({'$$\Delta \nu$$ (kHz)'},numPlots,1);
% 
% % subplot(3,1,3);
% ax = gca;
% 
% lineColors = {[120 137 163]/255,[33 54 86]/255,[133,125,245]/255,[0 0 0]/255};
% 
% hold on
% for ii = 1:numFilesIn
%     plot(dataOUT(:,1,ii),...
%         movmean((dataOUT(:,3,ii)-meanVal.freq(ii))/1e3,avg),...
%         'LineWidth',2,...
%         'Color',lineColors{ii+1});
% end
% hold off
% 
% xlim([min(dataOUT(:,1,1)) max(dataOUT(:,1,1))]);
% ylim([...
%     ((minVal.freq(2) - (span.freq(2)*.2/2)) - meanVal.freq(2))/1e3...
%     ((maxVal.freq(2) + (span.freq(2)*.2/2)) - meanVal.freq(2))/1e3...
%     ]);
% legend({'$$f_{OOL}$$','$$f_{CEP}$$'},'Location','northeast','interpreter','latex');
% 
% 
% % Make it look nice
% title(titles{1},'FontSize',30,'interpreter','latex');
% xlabel(xLabels{1},'FontSize',24,'interpreter','latex');
% ylabel(yLabels{1},'FontSize',24,'interpreter','latex');
% 
% ax.FontSize = 44;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(98);

titles = {'Frequency Drift, $$f_{OOL}$$ and $$f_{CEP}$$'};
xLabels = repmat({'Time (min)'},numPlots,1);
yLabels = {'$$\Delta \nu$$ (Hz)','$$\Delta \nu$$ (kHz)'};
legLabels = {'$$f_{OOL}$$','$$f_{CEP}$$'};
yScale = [1e1,1e3];

for ii = 1:numFilesIn
    subplot(2,1,ii);
    ax = gca;
    
    lineColors = {[33 54 86]/255};
    
    plot(dataOUT(:,1,ii),...
        movmean((dataOUT(:,3,ii)-meanVal.freq(ii))/yScale(ii),avg),...
        'LineWidth',2,...
        'Color',lineColors{1});
    
    xlim([min(dataOUT(:,1,ii)) max(dataOUT(:,1,ii))]);
    ylim([...
        ((minVal.freq(ii)-(span.freq(ii)*.2/2))-meanVal.freq(ii))/yScale(ii)...
        ((maxVal.freq(ii)+(span.freq(ii)*.2/2))-meanVal.freq(ii))/yScale(ii)...
        ]);
    legend(legLabels{ii},'Location','northeast','interpreter','latex');
    
    
    % Make it look nice
    ylabel(yLabels{ii},'FontSize',24,'interpreter','latex');
    
    if ii == 1
        title(titles{1},'FontSize',30,'interpreter','latex');
    end
    
    ax.FontSize = 44;
%     fillFig(0,-0.0)
    
    if ii < numFilesIn
        ax.XTick = [];
    end
    
end
